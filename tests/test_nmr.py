"""Secondary-structure comparison, including the metric that hid a bug.

Raw percent agreement reported a result for months that a majority-class
predictor beat. These tests pin the behaviour that exposed it.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from saxsnmr.nmr import (
    CSI_THRESHOLD,
    compare_to_dssp,
    csi_assignment,
    matthews_corrcoef,
    shift_dispersion,
)


def test_mcc_is_one_for_perfect_agreement():
    labels = np.array(list("HHHEEECCC"))
    np.testing.assert_allclose(matthews_corrcoef(labels, labels), 1.0)


def test_mcc_is_zero_for_a_constant_predictor():
    """The case raw agreement cannot see.

    A model answering "H" everywhere scores 50% agreement against this truth
    while containing no information. MCC reports 0, which is the point.
    """
    truth = np.array(list("HHHHEECC"))
    always_helix = np.array(["H"] * 8)
    agreement = float((truth == always_helix).mean())
    assert agreement == 0.5, "the trivial predictor should look respectable"
    mcc = matthews_corrcoef(truth, always_helix)
    assert abs(mcc) < 1e-9, f"MCC should be 0 for no information, got {mcc}"


def test_mcc_is_negative_when_predictions_are_anticorrelated():
    truth = np.array(list("HHHEEE"))
    flipped = np.array(list("EEEHHH"))
    assert matthews_corrcoef(truth, flipped) < -0.9


def _comparison(csi_labels, dssp_labels):
    resids = np.arange(1, len(csi_labels) + 1)
    table = pd.DataFrame({"resid": resids, "csi": list(csi_labels),
                          "d_diff": np.zeros(len(csi_labels))})
    dssp = np.array([list(dssp_labels)])          # one frame
    return compare_to_dssp(table, dssp, resids)


def test_comparison_reports_the_majority_baseline():
    """Every agreement figure must arrive with what a trivial model scores."""
    cmp = _comparison("HHHHHHEECC", "HHHHHHEECC")
    np.testing.assert_allclose(cmp.agreement, 1.0)
    # Six of ten residues are helix, so answering "H" always scores 60%.
    np.testing.assert_allclose(cmp.baseline_agreement, 0.6)
    assert cmp.beats_baseline


def test_a_result_below_its_own_baseline_is_flagged():
    """The failure that went unnoticed: respectable-looking agreement that a
    no-information predictor beats."""
    cmp = _comparison("HHHHHHHHEC", "CCCCCCCCEC")
    assert cmp.baseline_agreement == pytest.approx(0.8)
    assert cmp.agreement < cmp.baseline_agreement
    assert not cmp.beats_baseline
    assert "below majority baseline" in repr(cmp)


def test_confusion_matrix_shows_how_the_two_disagree():
    cmp = _comparison("HHHEEE", "HHHCCC")
    assert cmp.confusion.loc["H", "H"] == 3
    assert cmp.confusion.loc["E", "C"] == 3


def test_csi_threshold_is_on_the_difference_index_scale():
    """0.7 ppm is the textbook value for a single index; this code applies the
    threshold to dCA - dCB, which has about twice the dynamic range.

    Pinning this stops the earlier mistake from coming back: at 0.7 the coil
    band is a sliver and almost everything is forced into H or E.
    """
    assert 1.4 <= CSI_THRESHOLD <= 2.5, (
        "a threshold near 0.7 belongs to the single-index convention")


def test_csi_assignment_respects_its_threshold():
    sec = pd.DataFrame({
        "resid": [1, 2, 3, 4],
        "resname": ["ALA"] * 4,
        "d_ca": [0.0] * 4, "d_cb": [0.0] * 4,
        "d_diff": [5.0, -5.0, 0.0, np.nan],
    })
    labels = csi_assignment(sec)["csi"].tolist()
    assert labels == ["H", "E", "C", "-"]


def test_a_narrow_threshold_starves_the_coil_class():
    """Demonstrates the original bug rather than merely asserting the fix."""
    rng = np.random.default_rng(0)
    sec = pd.DataFrame({
        "resid": np.arange(200),
        "resname": ["ALA"] * 200,
        "d_ca": np.zeros(200), "d_cb": np.zeros(200),
        "d_diff": rng.normal(scale=4.6, size=200),   # measured spread
    })
    narrow = csi_assignment(sec, threshold=0.7)["csi"].value_counts()
    correct = csi_assignment(sec, threshold=1.75)["csi"].value_counts()
    assert narrow.get("C", 0) < correct.get("C", 0), (
        "the narrow threshold must assign fewer coil residues")
    assert narrow.get("C", 0) / 200 < 0.15, "coil band should be a sliver at 0.7"


def test_shift_dispersion_separates_informative_entries():
    """An entry with no conformational spread is not reporting structure.

    Selecting on this is safe because it looks only at the NMR data; it was
    added after an entry with 0.59 ppm spread produced an MCC of 0.009.
    """
    rng = np.random.default_rng(1)
    informative = pd.DataFrame({"d_diff": rng.normal(scale=4.0, size=120)})
    flat = pd.DataFrame({"d_diff": rng.normal(scale=0.5, size=120)})
    assert shift_dispersion(informative) > 2.0
    assert shift_dispersion(flat) < 2.0
