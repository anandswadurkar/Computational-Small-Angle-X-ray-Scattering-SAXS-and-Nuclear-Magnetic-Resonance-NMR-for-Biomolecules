"""Verification against closed-form answers.

These are the strongest tests in the suite: a solid sphere has an analytic
radius of gyration and an analytic scattering form factor, so agreement here
checks the Debye sum, the distance histogramming and the Guinier fit against
mathematics rather than against another implementation or deposited data.

Verification in the ASME sense -- is the code solving the equations correctly
-- as distinct from validation, which the notebooks do against experiment.
"""

from __future__ import annotations

import numpy as np
import pytest

from saxsnmr.saxs import ScatteringTables, guinier


def sphere_points(n: int, radius_A: float, seed: int = 0) -> np.ndarray:
    """`n` points uniform in a solid sphere of the given radius, in Angstrom.

    Rejection sampling rather than a radial transform, because getting the
    r^2 weighting wrong would bias Rg and quietly weaken every test here.
    """
    rng = np.random.default_rng(seed)
    out = np.empty((0, 3))
    while len(out) < n:
        batch = rng.uniform(-radius_A, radius_A, size=(4 * n, 3))
        inside = batch[(batch**2).sum(axis=1) <= radius_A**2]
        out = np.vstack([out, inside])
    return out[:n]


def sphere_form_factor(q: np.ndarray, radius_A: float) -> np.ndarray:
    """|F(q)|^2 for a uniform sphere, normalised to 1 at q = 0.

    F(q) = 3 [sin(qR) - qR cos(qR)] / (qR)^3
    """
    qr = np.asarray(q, dtype=float) * radius_A
    f = np.ones_like(qr)
    nz = qr > 1e-9
    f[nz] = 3.0 * (np.sin(qr[nz]) - qr[nz] * np.cos(qr[nz])) / qr[nz] ** 3
    return f**2


def uniform_tables(xyz: np.ndarray, bin_width: float = 0.5) -> ScatteringTables:
    """Histograms with every point treated as one identical scatterer.

    Uses a single-element alphabet so the profile below is the bare Debye
    sum with unit weights -- no form factors, no solvent term.
    """
    n = len(xyz)
    return ScatteringTables.from_coordinates(
        xyz, ["C"] * n, np.zeros(n), bin_width=bin_width, alphabet=("C",)
    )


def debye_unit_weights(tables: ScatteringTables, q: np.ndarray) -> np.ndarray:
    """I(q) = sum_ij sinc(q r_ij), the Debye sum with all weights 1."""
    return tables.sinc_matrix(q) @ tables.pair[("C", "C")]


RADIUS_A = 40.0
# Rg of a uniform solid sphere is R * sqrt(3/5); this is exact.
RG_EXACT = RADIUS_A * np.sqrt(3.0 / 5.0)


@pytest.fixture(scope="module")
def sphere():
    xyz = sphere_points(3000, RADIUS_A)
    return xyz, uniform_tables(xyz)


def test_sphere_radius_of_gyration_matches_closed_form(sphere):
    """The sampled point cloud reproduces R*sqrt(3/5) to sampling accuracy."""
    xyz, _ = sphere
    rg = np.sqrt(((xyz - xyz.mean(axis=0)) ** 2).sum(axis=1).mean())
    # 3000 points gives roughly 1% on Rg; 2% is a comfortable bound that
    # still fails loudly if the sampling becomes non-uniform.
    np.testing.assert_allclose(rg, RG_EXACT, rtol=0.02)


def test_debye_sum_reproduces_the_sphere_form_factor(sphere):
    """The histogrammed Debye sum matches |F(q)|^2 through the first minimum.

    Tested out to qR ~ 4.5, past the first form-factor minimum near qR = 4.49,
    so this checks the shape and not merely the q -> 0 limit.
    """
    _, tables = sphere
    q = np.linspace(0.005, 4.5 / RADIUS_A, 60)
    model = debye_unit_weights(tables, q)
    model = model / model[0]
    exact = sphere_form_factor(q, RADIUS_A)
    exact = exact / exact[0]
    # Finite sampling and 0.5 A binning put the floor here; the curve spans
    # two orders of magnitude, so a 6% envelope is a real constraint.
    np.testing.assert_allclose(model, exact, rtol=0.06, atol=0.01)


def test_guinier_recovers_the_analytic_rg(sphere):
    """Fitting the computed curve returns R*sqrt(3/5) without being told it."""
    _, tables = sphere
    q = np.linspace(0.002, 1.1 / RG_EXACT, 120)
    fit = guinier(q, debye_unit_weights(tables, q), None, min_points=25)
    assert fit.valid, f"Guinier fit rejected its own range: {fit.note}"
    np.testing.assert_allclose(fit.rg, RG_EXACT, rtol=0.03)
    assert fit.q_rg_max <= 1.3 + 1e-9, "fit exceeded its own validity limit"


def test_guinier_flags_a_particle_too_large_for_the_q_range():
    """A decaying curve whose Rg puts every point past qRg = 1.3 is flagged.

    This is the aggregate/fibril case: the fit succeeds arithmetically but
    outside the approximation's validity, so a number comes back with
    `valid` False rather than being quoted as an Rg.
    """
    q = np.linspace(0.01, 0.05, 40)
    huge_rg = 300.0  # qRg is already 3 at the first point
    decaying = np.exp(-(huge_rg**2) * q**2 / 3.0)
    fit = guinier(q, decaying, None, min_points=8)
    assert not fit.valid
    assert "no window" in fit.note
    assert fit.q_rg_min > 1.3, "the flag fired for the wrong reason"


def test_guinier_raises_when_the_curve_never_decays():
    """A curve rising with q has no radius of gyration at all.

    Distinct from the case above: there, an Rg exists but falls outside the
    valid range; here there is nothing to report, so returning a flagged
    number would invent one.
    """
    q = np.linspace(0.01, 0.2, 40)
    with pytest.raises(ValueError, match="does not decay"):
        guinier(q, np.exp(5.0 * q), None, min_points=8)


def test_scattering_tables_refuse_unknown_elements():
    """Silently dropping atoms would give a plausible curve for the wrong
    molecule, which is worse than failing."""
    xyz = sphere_points(50, 10.0)
    with pytest.raises(ValueError, match="no scattering parameters"):
        ScatteringTables.from_coordinates(xyz, ["C"] * 49 + ["Fe"], np.zeros(50))


def test_averaging_identical_frames_is_a_no_op(sphere):
    """Ensemble averaging must not change a single repeated conformation.

    If this fails, an ensemble chi-squared is wrong for reasons that have
    nothing to do with the physics being modelled.
    """
    _, tables = sphere
    q = np.linspace(0.01, 0.3, 30)
    one = debye_unit_weights(tables, q)
    averaged = ScatteringTables.average([tables, tables, tables])
    np.testing.assert_allclose(debye_unit_weights(averaged, q), one, rtol=1e-12)


@pytest.mark.parametrize("bin_width", [0.25, 0.5, 1.0])
def test_rg_is_insensitive_to_bin_width(bin_width):
    """Binning is a numerical choice, so the answer must not depend on it.

    A convergence check in the discretisation parameter: if Rg moved with
    bin width, the histogram step would be introducing physics.
    """
    xyz = sphere_points(2000, RADIUS_A, seed=1)
    tables = uniform_tables(xyz, bin_width=bin_width)
    q = np.linspace(0.002, 1.1 / RG_EXACT, 100)
    fit = guinier(q, debye_unit_weights(tables, q), None, min_points=25)
    np.testing.assert_allclose(fit.rg, RG_EXACT, rtol=0.04)
