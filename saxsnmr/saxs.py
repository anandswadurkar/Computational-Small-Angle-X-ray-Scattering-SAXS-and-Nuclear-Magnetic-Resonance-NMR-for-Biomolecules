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

__all__ = ["ScatteringTables", "profile", "fit", "SaxsFit"]

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
                         r_max: float | None = None) -> "ScatteringTables":
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

        d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
        limit = d.max() + bin_width if r_max is None else r_max + bin_width
        if d.max() > limit:
            raise ValueError(
                f"r_max={r_max} is smaller than the largest interatomic distance "
                f"({d.max():.1f} A); frames would bin inconsistently")
        edges = np.arange(0.0, limit, bin_width)
        centres = 0.5 * (edges[:-1] + edges[1:])
        flat = d.ravel()

        idx = {e: np.flatnonzero(np.array(elements) == e) for e in _ELEMENTS}
        idx = {e: v for e, v in idx.items() if v.size}

        pair: dict[tuple[str, str], np.ndarray] = {}
        for i, ea in enumerate(idx):
            for eb in list(idx)[i:]:
                sub = d[np.ix_(idx[ea], idx[eb])].ravel()
                h, _ = np.histogram(sub, bins=edges)
                # off-diagonal element pairs appear twice in the full double sum
                pair[(ea, eb)] = h.astype(np.float64) * (1.0 if ea == eb else 2.0)

        elem_sasa: dict[str, np.ndarray] = {}
        for ea, rows in idx.items():
            w = np.repeat(np.ones(rows.size), n) * np.tile(sasa_frac, rows.size)
            sub = d[rows, :].ravel()
            h, _ = np.histogram(sub, bins=edges, weights=w)
            elem_sasa[ea] = 2.0 * h            # cross term enters twice

        ww = np.outer(sasa_frac, sasa_frac).ravel()
        sasa_sasa, _ = np.histogram(flat, bins=edges, weights=ww)

        return cls(r=centres, pair=pair, elem_sasa=elem_sasa,
                   sasa_sasa=sasa_sasa.astype(np.float64), n_atoms=n)

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
                   q=q, model=scale * model, experiment=i_exp, sigma=sigma)
