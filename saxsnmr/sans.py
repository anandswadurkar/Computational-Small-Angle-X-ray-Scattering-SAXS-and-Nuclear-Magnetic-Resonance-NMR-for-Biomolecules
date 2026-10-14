"""Small-angle neutron scattering from coordinates.

SANS is computationally *easier* than SAXS here. Neutrons scatter from nuclei,
which are point scatterers, so the scattering length ``b`` carries no q
dependence -- there is no Cromer-Mann expansion to evaluate. The same Debye
sum applies,

    I(q) = sum_ij (b_i - rho_s V_i)(b_j - rho_s V_j) sinc(q r_ij)

and the distance histograms in ``ScatteringTables`` are reused unchanged.

What makes neutron scattering *physically* different is hydrogen:

* **Hydrogen scatters negatively** (b = -3.74 fm) where deuterium is strongly
  positive (+6.67 fm). Nothing in X-ray scattering behaves like this, and it
  is the whole basis of contrast variation.
* **Labile hydrogens exchange with the solvent.** Hydrogens on N, O and S
  swap with the buffer, so in D2O they become deuterium while the carbon-bound
  ones do not. The two classes therefore need different scattering lengths,
  which is why ``sans_tables`` builds histograms over an alphabet that
  distinguishes them.
* **Hydrogen scatters incoherently**, 80.3 barn against 1.76 barn coherent.
  That contributes a large flat background carrying no structural information,
  so a SANS fit must include a constant term. Omitting it forces the shape to
  absorb a baseline offset and biases Rg.

Because protein and solvent scattering-length densities cross near 40-45% D2O,
a measurement near that point has very little contrast and is dominated by the
hydration shell -- which this implementation does not model. The fitted
background and scale will happily absorb a great deal there, so results near
the match point should be treated as unreliable rather than merely noisy.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares, nnls

from .saxs import ScatteringTables, _DISPLACED_VOLUME

__all__ = ["NEUTRON_B", "solvent_sld", "exchangeable_hydrogens",
           "sans_elements", "sans_tables", "sans_profile", "sans_fit",
           "SansFit", "SANS_ALPHABET", "match_point_d2o",
           "scattering_weighted_rg"]

# Bound coherent scattering lengths, femtometres (Sears, Neutron News 1992).
# "Hx" is an exchangeable hydrogen: its value is set by the D2O fraction and
# is filled in at profile time, not here.
NEUTRON_B = {
    "H": -3.7390,
    "D": 6.671,
    "C": 6.6460,
    "N": 9.360,
    "O": 5.803,
    "P": 5.13,
    "S": 2.847,
}

# Incoherent cross sections, barn. Only hydrogen matters; everything else is
# negligible by comparison, which is why protonated buffers are so noisy.
INCOHERENT_BARN = {"H": 80.27, "D": 2.05, "C": 0.001, "N": 0.50,
                   "O": 0.0008, "P": 0.005, "S": 0.007}

SANS_ALPHABET = ("H", "Hx", "C", "N", "O", "P", "S")

V_WATER = 30.0          # A^3 per water molecule, as in the X-ray model

# Scattering length densities of the two solvents, fm per A^3.
SLD_H2O = (2 * NEUTRON_B["H"] + NEUTRON_B["O"]) / V_WATER
SLD_D2O = (2 * NEUTRON_B["D"] + NEUTRON_B["O"]) / V_WATER


def solvent_sld(d2o_fraction: float) -> float:
    """Scattering length density of an H2O/D2O mixture, fm/A^3.

    Linear in the D2O fraction, which is exact for a mixture of the two.
    """
    if not 0.0 <= d2o_fraction <= 1.0:
        raise ValueError(f"d2o_fraction must be in [0, 1], got {d2o_fraction}")
    return (1.0 - d2o_fraction) * SLD_H2O + d2o_fraction * SLD_D2O


def _b_exchangeable(d2o_fraction: float) -> float:
    """Mean scattering length of a labile hydrogen in a given solvent.

    Assumes labile positions equilibrate with the solvent, so the population
    is the solvent's H/D ratio. Buried amides exchange slowly and some never
    do, so this overstates deuteration for a folded protein -- the usual
    correction is an empirical exchange factor, which is not applied here.
    """
    return ((1.0 - d2o_fraction) * NEUTRON_B["H"]
            + d2o_fraction * NEUTRON_B["D"])


def exchangeable_hydrogens(topology) -> np.ndarray:
    """Boolean mask of hydrogens bonded to N, O or S.

    Those positions exchange with solvent; carbon-bound hydrogens do not.
    Determined from the bond graph rather than from atom names, because
    naming conventions vary between structures while bonding does not.
    """
    n_atoms = topology.n_atoms
    mask = np.zeros(n_atoms, dtype=bool)
    heavy_partner = {}
    for a, b in topology.bonds:
        for h, other in ((a, b), (b, a)):
            if h.element is not None and h.element.symbol == "H":
                heavy_partner[h.index] = other
    for h_index, other in heavy_partner.items():
        sym = other.element.symbol if other.element is not None else ""
        if sym in ("N", "O", "S"):
            mask[h_index] = True
    return mask


def sans_elements(topology) -> list[str]:
    """Per-atom labels for the neutron alphabet, marking labile H as ``Hx``."""
    labile = exchangeable_hydrogens(topology)
    out = []
    for atom in topology.atoms:
        sym = atom.element.symbol if atom.element is not None else ""
        if sym == "H":
            out.append("Hx" if labile[atom.index] else "H")
        else:
            out.append(sym)
    return out


def sans_tables(traj, frame: int = 0, bin_width: float = 0.25,
                r_max: float | None = None) -> ScatteringTables:
    """Distance histograms over the neutron alphabet for one frame.

    Separate from ``frame_tables`` because the alphabet differs: hydrogen is
    split into exchangeable and non-exchangeable, which the X-ray code has no
    reason to distinguish.

    SASA is still computed and stored, so a hydration term could be added
    later, but ``sans_profile`` does not currently use it.
    """
    import mdtraj as mdt

    xyz = traj.xyz[frame] * 10.0                     # nm -> Angstrom
    labels = sans_elements(traj.topology)
    sasa = mdt.shrake_rupley(traj[frame], mode="atom")[0]
    # Normalised the same way the X-ray path does, so the stored SASA tables
    # mean the same thing in both.
    from . import _MAX_SASA
    denom = np.array([_MAX_SASA.get(
        (a.element.symbol if a.element is not None else ""), 1.0)
        for a in traj.topology.atoms])
    frac = np.clip(sasa / denom, 0.0, 1.0)
    return ScatteringTables.from_coordinates(
        xyz, labels, frac, bin_width=bin_width, r_max=r_max,
        alphabet=SANS_ALPHABET)


def _weights(d2o_fraction: float, c1: float,
             protein_deuteration: float = 0.0) -> dict[str, float]:
    """Contrast-weighted scattering length per class, fm.

    ``protein_deuteration`` is the fraction of *non-exchangeable* hydrogen
    replaced by deuterium, i.e. biosynthetic labelling. This is not optional
    for real data: most deposited biomolecular SANS uses deuterated protein,
    because labelling is the reason to choose neutrons at all. Treating a
    perdeuterated sample as protonated gets the sign of the hydrogen
    contribution wrong for roughly a third of the atoms.

    Exchangeable positions follow the solvent regardless of labelling, since
    they re-equilibrate with the buffer.
    """
    if not 0.0 <= protein_deuteration <= 1.0:
        raise ValueError("protein_deuteration must be in [0, 1], got "
                         f"{protein_deuteration}")
    rho = solvent_sld(d2o_fraction)
    b_hx = _b_exchangeable(d2o_fraction)
    b_h = ((1.0 - protein_deuteration) * NEUTRON_B["H"]
           + protein_deuteration * NEUTRON_B["D"])
    out = {}
    for label in SANS_ALPHABET:
        element = "H" if label == "Hx" else label
        if label == "Hx":
            b = b_hx
        elif label == "H":
            b = b_h
        else:
            b = NEUTRON_B[element]
        out[label] = b - c1 * rho * _DISPLACED_VOLUME[element]
    return out


def sans_profile(tables: ScatteringTables, q: np.ndarray,
                 d2o_fraction: float, c1: float = 1.0,
                 sinc: np.ndarray | None = None,
                 protein_deuteration: float = 0.0) -> np.ndarray:
    """I(q) in arbitrary units, excluding the incoherent background.

    ``c1`` scales the displaced solvent volume exactly as in the X-ray model.
    No hydration shell term is included -- see the module docstring for why
    that matters most near the contrast match point.
    """
    q = np.asarray(q, dtype=np.float64)
    if sinc is None:
        sinc = tables.sinc_matrix(q)
    w = _weights(d2o_fraction, c1, protein_deuteration)
    intensity = np.zeros_like(q)
    for (ea, eb), hist in tables.pair.items():
        intensity += w[ea] * w[eb] * (sinc @ hist)
    return intensity


def match_point_d2o(tables: ScatteringTables,
                    protein_deuteration: float = 0.0) -> float:
    """D2O fraction at which this molecule's mean contrast vanishes.

    Solves <b> = rho_s(x) * <V> for x. Near this point the coherent signal is
    weakest and the hydration shell dominates what is left, so it marks where
    this implementation is least trustworthy rather than a target to aim for.
    """
    counts = {e: float(h.sum()) for e, h in tables.elem_sasa.items()}
    total = sum(counts.values())
    if total <= 0:
        return float("nan")

    def net(x: float) -> float:
        w = _weights(x, 1.0, protein_deuteration)
        return sum(n * w[label] for label, n in counts.items()) / total

    lo, hi = 0.0, 1.0
    if net(lo) * net(hi) > 0:
        return float("nan")          # no sign change: no match point
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if net(lo) * net(mid) <= 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


@dataclass
class SansFit:
    chi2: float
    scale: float
    background: float          # flat incoherent term, >= 0 by construction
    c1: float
    d2o_fraction: float
    protein_deuteration: float
    q: np.ndarray
    model: np.ndarray
    experiment: np.ndarray
    sigma: np.ndarray
    tables: "ScatteringTables | None" = None

    @property
    def background_fraction(self) -> float:
        """How much of the mean measured intensity the flat term accounts for.

        A large value means the fit is dominated by a featureless baseline and
        the structural part is carrying little of the signal.
        """
        mean = float(np.mean(self.experiment))
        return float(self.background / mean) if mean else float("nan")

    def __repr__(self) -> str:
        return (f"SansFit(chi2={self.chi2:.3f}, scale={self.scale:.4g}, "
                f"bg={self.background:.4g} ({self.background_fraction:.1%} of "
                f"mean I), c1={self.c1:.3f}, D2O={self.d2o_fraction:.0%}, "
                f"dProt={self.protein_deuteration:.0%})")


def sans_fit(tables: ScatteringTables, q: np.ndarray, i_exp: np.ndarray,
             sigma: np.ndarray | None, d2o_fraction: float,
             fit_background: bool = True,
             protein_deuteration: float = 0.0) -> SansFit:
    """Fit scale, a flat background and c1 to an experimental neutron curve.

    The background is not optional in practice. Incoherent scattering from
    hydrogen is large and flat, and a fit without it forces the structural
    shape to absorb a baseline offset, which biases Rg low. ``fit_background``
    exists so that effect can be demonstrated rather than assumed.

    Three free parameters against a SANS curve -- which typically has far
    fewer points and larger errors than a SAXS curve -- flatters a model even
    more than the X-ray case does. Report them alongside chi-squared, and read
    ``background_fraction`` before believing the shape.
    """
    q = np.asarray(q, dtype=np.float64)
    i_exp = np.asarray(i_exp, dtype=np.float64)
    if sigma is None:
        sigma = np.full_like(i_exp, np.abs(i_exp).mean() * 0.05)
    else:
        sigma = np.asarray(sigma, dtype=np.float64)
        bad = ~np.isfinite(sigma) | (sigma <= 0)
        if bad.any():                       # repair pointwise, as `fit` does
            good = np.abs(i_exp[~bad]).mean() * 0.05 if (~bad).any() else 1.0
            sigma = np.where(bad, good, sigma)

    sinc = tables.sinc_matrix(q)

    def model_for(c1: float) -> np.ndarray:
        return sans_profile(tables, q, d2o_fraction, c1=c1, sinc=sinc,
                            protein_deuteration=protein_deuteration)

    def residual(p):
        c1 = p[0]
        shape = model_for(c1)
        # Scale and background enter linearly, so solve for them exactly
        # rather than letting the optimiser wander over three dimensions.
        #
        # Both are constrained non-negative, which is physics rather than
        # convenience: incoherent scattering adds intensity and cannot
        # subtract it, and a negative scale would invert the model. An
        # unconstrained solve happily returns a negative background and then
        # the term is a nuisance parameter absorbing misfit, not a model of
        # anything.
        cols = [shape / sigma]
        if fit_background:
            cols.append(np.ones_like(q) / sigma)
        a = np.column_stack(cols)
        coef, _ = nnls(a, i_exp / sigma)
        resid = a @ coef - i_exp / sigma
        residual.coef = coef
        return resid

    sol = least_squares(residual, [1.0], bounds=([0.95], [1.05]),
                        xtol=1e-10, ftol=1e-10)
    c1 = float(sol.x[0])
    resid = residual([c1])
    coef = residual.coef
    scale = float(coef[0])
    background = float(coef[1]) if fit_background else 0.0
    n_par = 2 + (1 if fit_background else 0)
    dof = max(q.size - n_par, 1)
    chi2 = float((resid ** 2).sum() / dof)
    model = scale * model_for(c1) + background
    return SansFit(chi2=chi2, scale=scale, background=background, c1=c1,
                   d2o_fraction=d2o_fraction,
                   protein_deuteration=protein_deuteration, q=q, model=model,
                   experiment=i_exp, sigma=sigma, tables=tables)


def scattering_weighted_rg(traj, frame: int, d2o_fraction: float,
                           c1: float = 1.0,
                           protein_deuteration: float = 0.0) -> float:
    """Rg weighted by contrast-weighted scattering length, in Angstrom.

    **This is not the X-ray Rg.** Rg is a weighted second moment, and the
    weightings differ fundamentally: X-rays weight by electron count, where
    hydrogen is ~1 and everything is positive, while neutrons weight
    carbon-bound hydrogen *negatively*. Those hydrogens sit largely on
    outward-pointing side chains, so their negative contribution at large
    radius lowers the neutron-weighted Rg. For lysozyme in D2O the two differ
    by over an Angstrom, and that is physics rather than an error in either.

    This exists because it is the analytic value the Guinier limit of
    ``sans_profile`` must reproduce. Comparing the two validates the Debye
    sum, the scattering lengths and the contrast term together, against an
    answer computed a completely different way.

    Weights can be negative and may sum to near zero at the match point, where
    the weighted centroid and hence Rg stop being meaningful; NaN is returned
    rather than a large meaningless number.
    """
    xyz = np.asarray(traj.xyz[frame], dtype=np.float64) * 10.0
    labels = sans_elements(traj.topology)
    w_map = _weights(d2o_fraction, c1, protein_deuteration)
    w = np.array([w_map[l] for l in labels], dtype=np.float64)
    total = w.sum()
    if abs(total) < 1e-9 * np.abs(w).sum():
        return float("nan")
    centre = (w[:, None] * xyz).sum(axis=0) / total
    d2 = ((xyz - centre) ** 2).sum(axis=1)
    rg2 = float((w * d2).sum() / total)
    return float(np.sqrt(rg2)) if rg2 > 0 else float("nan")
