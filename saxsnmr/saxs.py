"""Compute SAXS profiles from coordinates and fit them to experimental curves.

Self-contained implementation of the standard three-term scattering model used
by CRYSOL and FoXS, written in numpy so the project has no binary or licensed
dependency. Each atom scatters as

    f_i(q) = f_atom(q) - c1 * f_displaced(q) + c2 * s_i * f_water(q)
             \\__vacuum__/   \\__bulk solvent__/   \\__hydration shell__/

with ``s_i`` the atom's fractional solvent accessibility, and the intensity
given by the Debye sum

    I(q) = sum_ij f_i(q) f_j(q) sinc(q r_ij)

Evaluating that sum directly costs O(N^2) per q value, which is far too slow to
repeat over a trajectory. Instead the q-independent geometry is precomputed once
as weighted distance histograms (see ``ScatteringTables``); afterwards a full
profile costs a single matrix product, so fitting c1/c2 and averaging over
frames are both cheap.

Fitting two free parameters will flatter almost any model, so ``fit`` returns
the fitted c1 and c2 alongside chi-squared. Report them together.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares

__all__ = ["ScatteringTables", "profile", "fit", "SaxsFit",
           "guinier", "GuinierFit", "pair_distribution", "rg_from_pr",
           "dmax_from_pr"]

# Cromer-Mann coefficients: f(q) = sum_k a_k exp(-b_k (q/4pi)^2) + c
_CROMER_MANN = {
    "H": ([0.493002, 0.322912, 0.140191, 0.040810],
          [10.5109, 26.1257, 3.14236, 57.7997], 0.003038),
    "C": ([2.310000, 1.020000, 1.588600, 0.865000],
          [20.8439, 10.2075, 0.568700, 51.6512], 0.215600),
    "N": ([12.21260, 3.132200, 2.012500, 1.166300],
          [0.005700, 9.893300, 28.99750, 0.582600], -11.52900),
    "O": ([3.048500, 2.286800, 1.546300, 0.867000],
          [13.27710, 5.701100, 0.323900, 32.90890], 0.250800),
    "P": ([6.434500, 4.179100, 1.780000, 1.490800],
          [1.906700, 27.15700, 0.526000, 68.16450], 1.114900),
    "S": ([6.905300, 5.203400, 1.437900, 1.586300],
          [1.467900, 22.21510, 0.253600, 56.17200], 0.866900),
}

# Solvent volume displaced per atom, cubic angstrom (Svergun et al.)
_DISPLACED_VOLUME = {"H": 5.15, "C": 16.44, "N": 2.49, "O": 9.13, "S": 19.86, "P": 5.73}

RHO_SOLVENT = 0.334       # bulk water electron density, e/A^3
V_WATER = 30.0            # volume of one water molecule, A^3
_ELEMENTS = tuple(_CROMER_MANN)


def _cdist(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Pairwise distances between two coordinate blocks.

    Uses scipy when available because it writes straight into the output array
    instead of materialising an (len(a), len(b), 3) difference first, which is
    the allocation this module is trying to avoid.
    """
    try:
        from scipy.spatial.distance import cdist
        return cdist(a, b)
    except Exception:
        return np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))


def _max_distance(xyz: np.ndarray, chunk: int) -> float:
    """Largest interatomic distance, computed without an N x N matrix."""
    best = 0.0
    for start in range(0, len(xyz), chunk):
        block = _cdist(xyz[start:start + chunk], xyz)
        best = max(best, float(block.max()))
    return best


def atomic_form_factor(element: str, q: np.ndarray) -> np.ndarray:
    a, b, c = _CROMER_MANN[element]
    s2 = (q / (4.0 * np.pi)) ** 2
    return sum(ai * np.exp(-bi * s2) for ai, bi in zip(a, b)) + c


def displaced_form_factor(element: str, q: np.ndarray) -> np.ndarray:
    """Gaussian sphere of bulk solvent excluded by the atom."""
    v = _DISPLACED_VOLUME[element]
    return RHO_SOLVENT * v * np.exp(-(v ** (2.0 / 3.0)) * q**2 / (4.0 * np.pi))


def water_form_factor(q: np.ndarray) -> np.ndarray:
    return RHO_SOLVENT * V_WATER * np.exp(-(V_WATER ** (2.0 / 3.0)) * q**2 / (4.0 * np.pi))


@dataclass
class ScatteringTables:
    """Distance histograms that depend only on geometry, not on q, c1 or c2.

    Expanding the Debye sum for f_i = u_Ei + c2 * s_i * w gives three groups of
    terms: element-pair, element-by-accessibility, and accessibility-squared.
    Each is a histogram over interatomic distance, computed once here.
    """

    r: np.ndarray                      # bin centres
    pair: dict[tuple[str, str], np.ndarray]
    elem_sasa: dict[str, np.ndarray]
    sasa_sasa: np.ndarray
    n_atoms: int

    @classmethod
    def from_coordinates(cls, xyz: np.ndarray, elements: list[str],
                         sasa_frac: np.ndarray, bin_width: float = 0.25,
                         r_max: float | None = None,
                         chunk: int = 256) -> "ScatteringTables":
        """Build tables for one conformation.

        ``r_max`` fixes the binning so tables from different frames can be
        averaged. Leave it unset for a single structure; set it to a value
        comfortably above the largest interatomic distance when processing a
        trajectory.
        """
        xyz = np.asarray(xyz, dtype=np.float64)
        sasa_frac = np.asarray(sasa_frac, dtype=np.float64)
        n = len(elements)
        if xyz.shape != (n, 3) or sasa_frac.shape != (n,):
            raise ValueError("coordinates, elements and sasa_frac disagree on atom count")

        # Distances are accumulated in row blocks rather than as one N x N
        # matrix. A 2000-atom protein would otherwise need a 92 MB (N, N, 3)
        # intermediate plus several 30 MB copies, and on a machine that is
        # short of memory that allocation fails intermittently -- killing the
        # interpreter rather than raising MemoryError. Blocking keeps the peak
        # near 5 MB and costs nothing in accuracy; the histograms are identical.
        limit = _max_distance(xyz, chunk) + bin_width if r_max is None else r_max + bin_width
        edges = np.arange(0.0, limit, bin_width)
        centres = 0.5 * (edges[:-1] + edges[1:])
        nbins = centres.size

        elements_arr = np.asarray(elements)
        idx = {e: np.flatnonzero(elements_arr == e) for e in _ELEMENTS}
        idx = {e: v for e, v in idx.items() if v.size}
        # position of each atom within its element group, for fast block slicing
        group_of = {e: np.zeros(n, dtype=bool) for e in idx}
        for e, rows in idx.items():
            group_of[e][rows] = True

        pair = {(ea, eb): np.zeros(nbins) for i, ea in enumerate(idx)
                for eb in list(idx)[i:]}
        elem_sasa = {e: np.zeros(nbins) for e in idx}
        sasa_sasa = np.zeros(nbins)
        order = list(idx)

        for start in range(0, n, chunk):
            stop = min(start + chunk, n)
            block = _cdist(xyz[start:stop], xyz)        # (chunk, n), ~4 MB
            if r_max is not None and block.max() > limit:
                raise ValueError(
                    f"r_max={r_max} is smaller than the largest interatomic distance "
                    f"({block.max():.1f} A); frames would bin inconsistently")
            rows_elem = elements_arr[start:stop]

            for ea in order:
                sel = np.flatnonzero(rows_elem == ea)
                if not sel.size:
                    continue
                sub = block[sel]                         # rows of element ea

                for eb in order:
                    key = (ea, eb) if (ea, eb) in pair else (eb, ea)
                    cols = idx[eb]
                    h, _ = np.histogram(sub[:, cols].ravel(), bins=edges)
                    # each unordered element pair is visited from both sides
                    # across the full sweep, which supplies the factor of two
                    # that off-diagonal pairs need; same-element pairs are
                    # visited once per ordered pair, which is also correct.
                    pair[key] += h

                w = np.broadcast_to(sasa_frac, sub.shape)
                h, _ = np.histogram(sub.ravel(), bins=edges, weights=w.ravel())
                elem_sasa[ea] += 2.0 * h

            ws = sasa_frac[start:stop, None] * sasa_frac[None, :]
            h, _ = np.histogram(block.ravel(), bins=edges, weights=ws.ravel())
            sasa_sasa += h
            del block

        return cls(r=centres, pair=pair, elem_sasa=elem_sasa,
                   sasa_sasa=sasa_sasa, n_atoms=n)

    @classmethod
    def average(cls, tables: list["ScatteringTables"]) -> "ScatteringTables":
        """Ensemble-average a set of tables built on identical binning.

        At fixed c1 and c2 the intensity is linear in these histograms, so the
        average of the histograms gives exactly the ensemble-averaged I(q).
        That is the quantity a scattering experiment measures, so the averaging
        must happen here -- before fitting -- not by averaging per-frame fits.
        """
        if not tables:
            raise ValueError("no tables to average")
        ref = tables[0]
        if any(t.r.shape != ref.r.shape for t in tables):
            raise ValueError("tables were built with different binning; pass a common r_max")
        n = float(len(tables))

        keys = set().union(*(set(t.pair) for t in tables))
        pair = {k: sum(t.pair.get(k, 0.0) for t in tables) / n for k in keys}
        ekeys = set().union(*(set(t.elem_sasa) for t in tables))
        elem = {k: sum(t.elem_sasa.get(k, 0.0) for t in tables) / n for k in ekeys}
        ss = sum(t.sasa_sasa for t in tables) / n
        return cls(r=ref.r, pair=pair, elem_sasa=elem, sasa_sasa=ss, n_atoms=ref.n_atoms)

    def sinc_matrix(self, q: np.ndarray) -> np.ndarray:
        """sinc(q r) for every (q, bin) pair; the only O(nq * nbins) step."""
        qr = np.outer(q, self.r)
        out = np.ones_like(qr)
        nz = qr > 1e-12
        out[nz] = np.sin(qr[nz]) / qr[nz]
        return out


def profile(tables: ScatteringTables, q: np.ndarray, c1: float = 1.0,
            c2: float = 0.0, sinc: np.ndarray | None = None) -> np.ndarray:
    """I(q) for the given excluded-volume and hydration parameters."""
    q = np.asarray(q, dtype=np.float64)
    if sinc is None:
        sinc = tables.sinc_matrix(q)

    u = {e: atomic_form_factor(e, q) - c1 * displaced_form_factor(e, q)
         for e in tables.elem_sasa}
    w = water_form_factor(q)

    intensity = np.zeros_like(q)
    for (ea, eb), hist in tables.pair.items():
        intensity += u[ea] * u[eb] * (sinc @ hist)
    if c2:
        for ea, hist in tables.elem_sasa.items():
            intensity += c2 * u[ea] * w * (sinc @ hist)
        intensity += (c2 * w) ** 2 * (sinc @ tables.sasa_sasa)
    return intensity


# Electrons per atom, for weighting the real-space distance distribution.
_Z = {"H": 1, "C": 6, "N": 7, "O": 8, "P": 15, "S": 16}


def pair_distribution(tables: ScatteringTables,
                      normalise: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Electron-weighted distance distribution P(r) from precomputed tables.

    P(r) is the q=0 real-space view of the same geometry the Debye sum uses:
    every interatomic distance, weighted by the product of the two atoms'
    electron counts. It is returned straight from the element-pair histograms,
    so it costs nothing once the tables exist.

    This is the model's own P(r), *not* an indirect Fourier transform of the
    measured intensity. Comparing it to an experimental GNOM P(r) compares two
    quantities computed in genuinely different ways -- which is the point, but
    it also means small differences in shape are not automatically a flaw in
    the model: GNOM's result depends on the chosen Dmax and regularisation.

    The r=0 bin is dropped. It holds each atom's distance to itself, which is
    an artifact of how the histogram is accumulated rather than a real pair.
    """
    p = np.zeros_like(tables.r, dtype=np.float64)
    for (ea, eb), hist in tables.pair.items():
        p += _Z[ea] * _Z[eb] * hist
    p[tables.r < 1e-9] = 0.0
    if normalise and p.max() > 0:
        p = p / p.max()
    return tables.r, p


def rg_from_pr(r: np.ndarray, p: np.ndarray) -> float:
    """Radius of gyration from the second moment of P(r).

    Rg^2 = integral r^2 P(r) dr / (2 integral P(r) dr). This uses the whole
    curve rather than the low-q limit, so it is an independent check on the
    Guinier value and does not share its sensitivity to the fitted range.
    """
    p = np.asarray(p, dtype=np.float64)
    total = p.sum()
    if total <= 0:
        return float("nan")
    return float(np.sqrt((r ** 2 * p).sum() / (2.0 * total)))


def dmax_from_pr(r: np.ndarray, p: np.ndarray, floor: float = 0.01) -> float:
    """Largest distance at which P(r) is still above `floor` of its peak.

    A threshold is needed because a computed P(r) has no sharp end: a handful
    of atom pairs sit beyond the bulk of the molecule. GNOM's Dmax, by
    contrast, is a fitted parameter. The two are not the same quantity and
    exact agreement is not expected.
    """
    p = np.asarray(p, dtype=np.float64)
    if p.max() <= 0:
        return float("nan")
    above = np.flatnonzero(p >= floor * p.max())
    return float(r[above[-1]]) if above.size else float("nan")


@dataclass
class GuinierFit:
    """Result of a Guinier fit, with the information needed to judge it.

    `rg_stderr` is the **formal** error from the fit covariance, which assumes
    the deposited sigmas are exactly right and that the Guinier law holds over
    the whole fitted window. It is routinely far smaller than the real
    uncertainty: on lysozyme `SASDMJ2` this gives 0.03 A where SASBDB's own
    ATSAS analysis quotes 0.26 A. The difference is not a bug in either -- it
    is what "formal error" means. Do not quote `rg_stderr` as an uncertainty on
    Rg; use it to compare fits computed the same way, and treat the systematic
    choice of fitting range as the dominant term.
    """

    rg: float                  # Angstrom
    rg_stderr: float           # formal only -- see the class docstring
    i0: float
    i0_stderr: float
    first: int                 # index range used, into the curve passed in
    last: int
    q: np.ndarray              # the q values actually fitted
    ln_intensity: np.ndarray
    ln_model: np.ndarray
    q_rg_min: float
    q_rg_max: float
    reduced_chi2: float        # of the straight line itself; >> 1 means curvature
    valid: bool
    note: str

    @property
    def n_points(self) -> int:
        return self.last - self.first

    def __repr__(self) -> str:
        flag = "" if self.valid else "  [INVALID]"
        return (f"GuinierFit(Rg={self.rg:.2f} A, I0={self.i0:.4g}, "
                f"n={self.n_points}, qRg={self.q_rg_min:.2f}-{self.q_rg_max:.2f}, "
                f"chi2red={self.reduced_chi2:.2f}){flag}")


def _guinier_window(q2: np.ndarray, ln_i: np.ndarray, w: np.ndarray,
                    first: int, last: int):
    """Weighted straight-line fit of ln I against q^2 over one window."""
    x, y, wi = q2[first:last], ln_i[first:last], w[first:last]
    design = np.column_stack([np.ones_like(x), x])
    a = design.T @ (wi[:, None] * design)
    b = design.T @ (wi * y)
    try:
        beta = np.linalg.solve(a, b)
        cov = np.linalg.inv(a)
    except np.linalg.LinAlgError:
        return None
    intercept, slope = beta
    if not np.isfinite(slope) or slope >= 0:
        return None                      # upward curve: no real Rg
    rg = float(np.sqrt(-3.0 * slope))
    slope_err = float(np.sqrt(max(cov[1, 1], 0.0)))
    model = intercept + slope * x
    dof = max(len(x) - 2, 1)
    chi2red = float((wi * (y - model) ** 2).sum() / dof)
    return {
        "rg": rg,
        "rg_stderr": float(1.5 / rg * slope_err) if rg > 0 else float("nan"),
        "i0": float(np.exp(intercept)),
        "i0_stderr": float(np.exp(intercept) * np.sqrt(max(cov[0, 0], 0.0))),
        "ln_model": model,
        "reduced_chi2": chi2red,
    }


def guinier(q: np.ndarray, intensity: np.ndarray,
            sigma: np.ndarray | None = None, q_rg_limit: float = 1.3,
            min_points: int = 8, skip: int = 0) -> GuinierFit:
    """Fit ln I(q) = ln I0 - Rg^2 q^2 / 3 over a self-consistent low-q range.

    The Guinier approximation holds only while q*Rg is small, conventionally
    below about 1.3 for a globular particle. That criterion is circular -- the
    valid range depends on the Rg being measured -- so the range is chosen by
    extending the window one point at a time and stopping when the fitted Rg
    puts the last point past the limit.

    Reporting Rg without the range it came from is meaningless, so the fitted
    indices, the qRg span and a validity flag are all returned. `valid` being
    False means the result failed its own criterion; the number is still there
    to look at, but it should not be quoted.

    `skip` drops leading points. Data very near the beamstop is often
    unreliable, and SASBDB curves vary in how much of it survived processing.
    """
    q = np.asarray(q, dtype=np.float64)
    intensity = np.asarray(intensity, dtype=np.float64)

    usable = np.isfinite(q) & np.isfinite(intensity) & (intensity > 0)
    if sigma is not None:
        sigma = np.asarray(sigma, dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            w = np.where(np.isfinite(sigma) & (sigma > 0),
                         (intensity / sigma) ** 2, 1.0)
    else:
        w = np.ones_like(q)
    w = np.where(usable, w, 0.0)

    q2 = q ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        ln_i = np.where(usable, np.log(np.abs(intensity)), 0.0)

    start = int(skip)
    n = q.size
    if n - start < min_points:
        raise ValueError(f"only {n - start} usable points; need {min_points}")

    best, best_last = None, None
    for last in range(start + min_points, n + 1):
        got = _guinier_window(q2, ln_i, w, start, last)
        if got is None:
            continue
        if q[last - 1] * got["rg"] <= q_rg_limit:
            best, best_last = got, last
        elif best is not None:
            break                        # past the valid range; keep the last good one

    valid = best is not None
    note = ""
    if not valid:
        # Fall back to the smallest window so there is something to inspect.
        best = _guinier_window(q2, ln_i, w, start, start + min_points)
        best_last = start + min_points
        if best is None:
            raise ValueError("Guinier fit failed: intensity does not decay at low q")
        note = (f"no window satisfies qRg <= {q_rg_limit}; showing the first "
                f"{min_points} points. The curve may start above the Guinier "
                f"region, or the particle may be aggregated.")
    elif start > 0:
        note = f"skipped {start} leading point(s)"

    used = slice(start, best_last)
    return GuinierFit(
        rg=best["rg"], rg_stderr=best["rg_stderr"],
        i0=best["i0"], i0_stderr=best["i0_stderr"],
        first=start, last=best_last,
        q=q[used], ln_intensity=ln_i[used], ln_model=best["ln_model"],
        q_rg_min=float(q[start] * best["rg"]),
        q_rg_max=float(q[best_last - 1] * best["rg"]),
        reduced_chi2=best["reduced_chi2"],
        valid=valid, note=note,
    )


@dataclass
class SaxsFit:
    chi2: float
    scale: float
    c1: float
    c2: float
    q: np.ndarray
    model: np.ndarray
    experiment: np.ndarray
    sigma: np.ndarray
    # The geometry the fit was computed from, kept so real-space quantities
    # can be derived without rebuilding it. For an ensemble fit this is the
    # averaged tables, which is otherwise not recoverable from the caller.
    tables: "ScatteringTables | None" = None

    def __repr__(self) -> str:
        return (f"SaxsFit(chi2={self.chi2:.3f}, scale={self.scale:.4g}, "
                f"c1={self.c1:.3f}, c2={self.c2:.3f}, n={self.q.size})")


def fit(tables: ScatteringTables, q: np.ndarray, i_exp: np.ndarray,
        sigma: np.ndarray | None = None, fit_hydration: bool = True) -> SaxsFit:
    """Fit scale, c1 and c2 to an experimental curve; report reduced chi-squared.

    c1 is bounded to 0.95-1.05 and c2 to 0-4.0, the ranges used by CRYSOL and
    FoXS. Without bounds these absorb genuine structural error and the fit stops
    meaning anything.
    """
    q = np.asarray(q, dtype=np.float64)
    i_exp = np.asarray(i_exp, dtype=np.float64)

    # Repair bad uncertainties pointwise. Discarding the whole array because a
    # few points are unusable would silently turn chi-squared into an unweighted
    # residual, which is a different and much less meaningful quantity.
    fallback = float(np.abs(i_exp).mean()) * 0.01
    if sigma is None:
        sigma = np.full_like(i_exp, fallback)
    else:
        sigma = np.array(sigma, dtype=np.float64, copy=True)
        bad = ~np.isfinite(sigma) | (sigma <= 0)
        if bad.any():
            good = sigma[~bad]
            sigma[bad] = float(np.median(good)) if good.size else fallback

    sinc = tables.sinc_matrix(q)

    def residual(p: np.ndarray) -> np.ndarray:
        c1, c2 = p[0], (p[1] if fit_hydration else 0.0)
        model = profile(tables, q, c1=c1, c2=c2, sinc=sinc)
        scale = np.sum(model * i_exp / sigma**2) / np.sum(model**2 / sigma**2)
        return (scale * model - i_exp) / sigma

    p0 = [1.0, 0.5] if fit_hydration else [1.0]
    lo = [0.95, 0.0] if fit_hydration else [0.95]
    hi = [1.05, 4.0] if fit_hydration else [1.05]
    sol = least_squares(residual, p0, bounds=(lo, hi), xtol=1e-10, ftol=1e-10)

    c1 = sol.x[0]
    c2 = sol.x[1] if fit_hydration else 0.0
    model = profile(tables, q, c1=c1, c2=c2, sinc=sinc)
    scale = np.sum(model * i_exp / sigma**2) / np.sum(model**2 / sigma**2)
    dof = max(q.size - (3 if fit_hydration else 2), 1)
    chi2 = float(np.sum(((scale * model - i_exp) / sigma) ** 2) / dof)

    return SaxsFit(chi2=chi2, scale=float(scale), c1=float(c1), c2=float(c2),
                   q=q, model=scale * model, experiment=i_exp, sigma=sigma,
                   tables=tables)
