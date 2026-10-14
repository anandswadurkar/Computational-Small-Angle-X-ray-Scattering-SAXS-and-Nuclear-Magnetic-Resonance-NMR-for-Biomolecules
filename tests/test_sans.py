"""Neutron forward model, against tabulated and analytic answers.

None of these needs network or a structure from a database: the scattering
lengths are nuclear data, and the Guinier theorem holds for any point cloud.
"""

from __future__ import annotations

import numpy as np
import pytest

from saxsnmr.sans import (
    NEUTRON_B,
    SLD_D2O,
    SLD_H2O,
    _b_exchangeable,
    _weights,
    solvent_sld,
)


def test_solvent_densities_match_published_values():
    """H2O and D2O scattering length densities, in fm/A^3.

    Quoted in every SANS reference as -0.56e-6 and 6.38e-6 per A^2, which is
    -0.056 and 0.638 fm/A^3 for a 30 A^3 water. Getting the sign of H2O wrong
    is the classic error, and it would be invisible in a fitted scale factor.
    """
    np.testing.assert_allclose(SLD_H2O, -0.0558, atol=0.002)
    np.testing.assert_allclose(SLD_D2O, 0.6382, atol=0.002)
    assert SLD_H2O < 0 < SLD_D2O, "the sign reversal is the basis of contrast"


def test_hydrogen_scatters_negatively_and_deuterium_does_not():
    """The one fact that makes neutron scattering different from X-ray."""
    assert NEUTRON_B["H"] < 0
    assert NEUTRON_B["D"] > 0
    assert abs(NEUTRON_B["D"]) > abs(NEUTRON_B["H"])


@pytest.mark.parametrize("fraction", [0.0, 0.25, 0.5, 0.75, 1.0])
def test_solvent_sld_is_linear_in_the_d2o_fraction(fraction):
    """Exact for a mixture of two solvents, so this is an identity, not a fit."""
    expected = (1 - fraction) * SLD_H2O + fraction * SLD_D2O
    np.testing.assert_allclose(solvent_sld(fraction), expected, rtol=1e-12)


@pytest.mark.parametrize("bad", [-0.01, 1.01, 2.0])
def test_solvent_sld_refuses_impossible_fractions(bad):
    with pytest.raises(ValueError, match="d2o_fraction"):
        solvent_sld(bad)


def test_exchangeable_hydrogen_follows_the_solvent():
    """Labile positions equilibrate with the buffer, so b moves H -> D."""
    np.testing.assert_allclose(_b_exchangeable(0.0), NEUTRON_B["H"])
    np.testing.assert_allclose(_b_exchangeable(1.0), NEUTRON_B["D"])
    midpoint = _b_exchangeable(0.5)
    assert NEUTRON_B["H"] < midpoint < NEUTRON_B["D"]


def test_protein_deuteration_moves_only_the_fixed_hydrogens():
    """Biosynthetic labelling changes carbon-bound H; exchangeable H already
    follow the solvent and must not be double-counted."""
    plain = _weights(0.0, 1.0, protein_deuteration=0.0)
    labelled = _weights(0.0, 1.0, protein_deuteration=1.0)
    assert labelled["H"] > plain["H"], "labelling must raise fixed-H scattering"
    np.testing.assert_allclose(labelled["Hx"], plain["Hx"], rtol=1e-12)
    np.testing.assert_allclose(labelled["C"], plain["C"], rtol=1e-12)


@pytest.mark.parametrize("bad", [-0.1, 1.5])
def test_weights_refuse_impossible_deuteration(bad):
    with pytest.raises(ValueError, match="protein_deuteration"):
        _weights(0.0, 1.0, protein_deuteration=bad)


# ---------------------------------------------------------------------------
# The Guinier theorem: the low-q limit of the Debye sum must reproduce the
# analytic scattering-length-weighted Rg. The two share no code beyond the
# coordinates, so this checks the scattering lengths, the contrast term and
# the histogramming together.
# ---------------------------------------------------------------------------

def _toy_molecule(n: int = 400, seed: int = 0):
    """A point cloud with a realistic mix of the neutron atom classes."""
    rng = np.random.default_rng(seed)
    xyz = rng.normal(scale=9.0, size=(n, 3))
    labels = rng.choice(["H", "Hx", "C", "N", "O"], size=n,
                        p=[0.38, 0.14, 0.30, 0.09, 0.09]).tolist()
    return xyz, labels


def _weighted_rg(xyz, labels, d2o, dprot=0.0):
    w_map = _weights(d2o, 1.0, dprot)
    w = np.array([w_map[x] for x in labels])
    total = w.sum()
    centre = (w[:, None] * xyz).sum(axis=0) / total
    return float(np.sqrt((w * ((xyz - centre) ** 2).sum(axis=1)).sum() / total))


@pytest.mark.parametrize("d2o,dprot", [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0),
                                       (0.2, 0.5)])
def test_guinier_limit_matches_the_analytic_weighted_rg(d2o, dprot):
    from saxsnmr.sans import SANS_ALPHABET, sans_profile
    from saxsnmr.saxs import ScatteringTables, guinier

    xyz, labels = _toy_molecule()
    tables = ScatteringTables.from_coordinates(
        xyz, labels, np.zeros(len(labels)), bin_width=0.25,
        alphabet=SANS_ALPHABET)
    analytic = _weighted_rg(xyz, labels, d2o, dprot)
    q = np.linspace(0.002, 0.9 / analytic, 140)
    fit = guinier(q, sans_profile(tables, q, d2o, protein_deuteration=dprot),
                  None, min_points=30)
    np.testing.assert_allclose(fit.rg, analytic, rtol=0.03)


def test_contrast_vanishes_somewhere_between_the_two_solvents():
    """A protonated protein must have a match point; I(0) collapses there."""
    from saxsnmr.sans import SANS_ALPHABET, match_point_d2o, sans_profile
    from saxsnmr.saxs import ScatteringTables

    xyz, labels = _toy_molecule()
    tables = ScatteringTables.from_coordinates(
        xyz, labels, np.zeros(len(labels)), alphabet=SANS_ALPHABET)
    mp = match_point_d2o(tables)
    assert 0.2 < mp < 0.8, f"match point {mp:.2f} is not physical"

    q = np.array([0.01])
    at_match = abs(sans_profile(tables, q, mp)[0])
    extremes = min(abs(sans_profile(tables, q, x)[0]) for x in (0.0, 1.0))
    assert at_match < 0.02 * extremes


def test_perdeuterated_protein_has_no_match_point():
    """It scatters more strongly than pure D2O, so no mixture can match it.

    This is why deuterium labelling works: a labelled component stays visible
    whatever the solvent, while its unlabelled partner can be matched out.
    """
    from saxsnmr.sans import SANS_ALPHABET, match_point_d2o
    from saxsnmr.saxs import ScatteringTables

    xyz, labels = _toy_molecule()
    tables = ScatteringTables.from_coordinates(
        xyz, labels, np.zeros(len(labels)), alphabet=SANS_ALPHABET)
    assert np.isfinite(match_point_d2o(tables, protein_deuteration=0.0))
    assert not np.isfinite(match_point_d2o(tables, protein_deuteration=1.0))


def test_fitted_background_is_never_negative():
    """Incoherent scattering adds intensity; it cannot subtract it.

    Without the non-negative constraint the term silently becomes a nuisance
    parameter absorbing misfit, which is how it behaved before nnls.
    """
    from saxsnmr.sans import SANS_ALPHABET, sans_fit, sans_profile
    from saxsnmr.saxs import ScatteringTables

    xyz, labels = _toy_molecule()
    tables = ScatteringTables.from_coordinates(
        xyz, labels, np.zeros(len(labels)), alphabet=SANS_ALPHABET)
    q = np.linspace(0.01, 0.3, 80)
    truth = sans_profile(tables, q, 1.0)
    # Data sitting *below* the model is what tempts the fit to go negative.
    observed = 0.5 * truth / truth[0]
    fit = sans_fit(tables, q, observed, None, 1.0)
    assert fit.background >= 0.0
    assert fit.scale >= 0.0
