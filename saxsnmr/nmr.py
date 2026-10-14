"""Compare NMR chemical shifts against secondary structure from a trajectory.

Predicting raw chemical shifts from coordinates needs SPARTA+, SHIFTX2 or
UCBShift. None of those has a usable Windows build, so this module validates
against NMR through the route that is fully portable: secondary chemical
shifts.

The 13C-alpha and 13C-beta shifts of a residue depend strongly and predictably
on its backbone conformation. Subtracting tabulated random-coil values leaves
the *secondary* shift, whose sign reports helix versus strand -- the chemical
shift index of Wishart and Sykes. That gives a per-residue, experimentally
grounded secondary structure assignment which can be compared directly with
DSSP computed over the MD ensemble.

This is a weaker test than a direct shift comparison: it checks that the
simulated ensemble has the right secondary structure where NMR says it should,
not that every shift is reproduced. It is honest about what it measures and it
runs anywhere.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "ShiftComparison",
    "best_offset",
    "compare_to_dssp",
    "csi_assignment",
    "matthews_corrcoef",
    "secondary_shifts",
    "shift_dispersion",
    "threshold_scan",
]

# Random-coil 13C shifts, ppm (Wishart & Sykes). Glycine has no C-beta.
_RANDOM_COIL_CA = {
    "ALA": 52.5, "ARG": 56.0, "ASN": 52.5, "ASP": 54.0, "CYS": 58.2,
    "GLN": 55.7, "GLU": 56.4, "GLY": 45.1, "HIS": 55.0, "ILE": 61.1,
    "LEU": 55.1, "LYS": 56.5, "MET": 55.4, "PHE": 57.7, "PRO": 63.3,
    "SER": 58.3, "THR": 61.8, "TRP": 57.5, "TYR": 57.9, "VAL": 62.2,
}
_RANDOM_COIL_CB = {
    "ALA": 19.0, "ARG": 30.9, "ASN": 38.9, "ASP": 41.0, "CYS": 28.0,
    "GLN": 29.4, "GLU": 29.7, "HIS": 29.0, "ILE": 38.8, "LEU": 42.4,
    "LYS": 33.1, "MET": 32.9, "PHE": 39.6, "PRO": 32.1, "SER": 63.8,
    "THR": 69.8, "TRP": 29.6, "TYR": 38.8, "VAL": 32.9,
}

# Threshold on the difference index dCA - dCB, in ppm.
#
# The canonical secondary shifts are helix dCA +3.1 / dCB -0.4, and strand
# dCA -1.5 / dCB +2.2. The difference index therefore separates helix (+3.5)
# from strand (-3.7) with roughly twice the dynamic range of either shift
# alone, because helix raises CA *and* lowers CB. Assigning H or E when the
# index exceeds half the canonical magnitude puts the boundary at ~1.75 ppm.
#
# This was previously 0.7, which is the textbook threshold for a *single*
# index applied by mistake to the difference. With a measured standard
# deviation of 4.65 ppm, a +/-0.7 band is a sliver through the middle of the
# distribution: on lysozyme it assigned coil to 10 residues out of 129 where
# DSSP finds 63, forcing almost everything into H or E. Agreement then fell
# below what a majority-class predictor achieves, which is what made the
# problem visible.
#
# Do not tune this against a reference structure. Choosing it to maximise
# agreement is fitting to the answer; `threshold_scan` exists so the
# sensitivity can be reported instead of hidden.
CSI_THRESHOLD = 1.75


def secondary_shifts(shifts):
    """Per-residue secondary shifts from a BMRB assigned-shift table.

    Expects the DataFrame returned by ``Fetcher.bmrb_shifts`` (columns
    ``resid``, ``resname``, ``atom``, ``shift``). Returns one row per residue
    with ``d_ca``, ``d_cb`` and their difference, which is the most robust
    single indicator.
    """
    import pandas as pd

    df = shifts.copy()
    df["resname"] = df["resname"].str.upper().str.strip()
    df["atom"] = df["atom"].str.upper().str.strip()
    wanted = df[df["atom"].isin(["CA", "CB"])]

    rows = []
    for (resid, resname), grp in wanted.groupby(["resid", "resname"], dropna=True):
        if resname not in _RANDOM_COIL_CA:
            continue
        ca = grp.loc[grp["atom"] == "CA", "shift"]
        cb = grp.loc[grp["atom"] == "CB", "shift"]
        d_ca = float(ca.iloc[0]) - _RANDOM_COIL_CA[resname] if len(ca) else np.nan
        d_cb = (float(cb.iloc[0]) - _RANDOM_COIL_CB[resname]
                if len(cb) and resname in _RANDOM_COIL_CB else np.nan)
        rows.append({"resid": int(resid), "resname": resname,
                     "d_ca": d_ca, "d_cb": d_cb,
                     "d_diff": d_ca - d_cb if np.isfinite(d_ca) and np.isfinite(d_cb)
                     else d_ca})
    out = pd.DataFrame(rows).sort_values("resid").reset_index(drop=True)
    return out


def shift_dispersion(sec_shifts) -> float:
    """Spread of the difference index, in ppm -- is this entry informative?

    A folded protein with mixed secondary structure puts helix near +3.5 ppm
    on this index and strand near -3.7, so the population standard deviation
    should be around 3 ppm or more. An entry well below ~2 ppm is not
    reporting conformational dispersion: it may be a denatured or
    dynamically averaged sample, or shifts that already had a reference
    subtracted before deposition.

    Selecting entries on this is safe, because it looks only at the NMR data.
    Selecting them on agreement with a structure would be fitting to the
    answer. Among the RNase A entries, BMRB 16742 gives 0.59 ppm where five
    others give 3.6-5.0 -- picking it produced an MCC of 0.009, and the
    controls caught that only after the fact.
    """
    import numpy as _np
    return float(_np.nanstd(sec_shifts["d_diff"].to_numpy(dtype=float)))


def csi_assignment(sec_shifts, threshold: float = CSI_THRESHOLD):
    """Label each residue H (helix), E (strand) or C (coil) from its shifts.

    ``threshold`` is in ppm on the difference index; see ``CSI_THRESHOLD`` for
    where the default comes from and why it should not be tuned against a
    reference structure.
    """

    def label(row) -> str:
        v = row["d_diff"]
        if not np.isfinite(v):
            return "-"
        if v > threshold:
            return "H"
        if v < -threshold:
            return "E"
        return "C"

    out = sec_shifts.copy()
    out["csi"] = out.apply(label, axis=1)
    return out


def matthews_corrcoef(y_true: np.ndarray, y_pred: np.ndarray,
                      labels: tuple[str, ...] = ("H", "E", "C")) -> float:
    """Multiclass Matthews correlation coefficient.

    0 is chance, 1 is perfect, negative is anti-correlated. Unlike raw
    agreement it accounts for class imbalance, which matters here: secondary
    structure classes are far from equally populated, so a model that always
    answers "helix" can post a high agreement while containing no information.
    """
    y_true, y_pred = np.asarray(y_true), np.asarray(y_pred)
    c = np.array([[np.sum((y_true == a) & (y_pred == b)) for b in labels]
                  for a in labels], dtype=float)
    t, p, n = c.sum(axis=1), c.sum(axis=0), c.sum()
    num = c.trace() * n - (t * p).sum()
    den = np.sqrt((n ** 2 - (p ** 2).sum()) * (n ** 2 - (t ** 2).sum()))
    if den <= 0:
        # One of the two labellings is constant, so MCC is formally
        # undefined. Zero is the right summary and NaN is not: a predictor
        # that answers the same class everywhere carries no information, and
        # that is exactly the degenerate case this metric exists to expose.
        # Returning NaN would propagate into a report as a blank rather than
        # as the verdict.
        return 0.0
    return float(num / den)


@dataclass
class ShiftComparison:
    """Agreement between NMR-derived and simulated secondary structure.

    Read ``mcc`` first. ``agreement`` is the fraction of residues matching,
    which is intuitive and misleading: compare it with ``baseline_agreement``,
    what a predictor that always answers the commonest class would score. On
    lysozyme that baseline is around 50-60%, so an agreement in that range
    means nothing on its own.
    """

    n_compared: int
    agreement: float               # fraction of residues matching H/E/C
    baseline_agreement: float      # majority-class predictor, for comparison
    balanced_accuracy: float       # mean per-class recall
    mcc: float                     # Matthews correlation; 0 = chance
    helix_agreement: float
    strand_agreement: float
    coil_agreement: float
    confusion: object              # DataFrame, NMR class x simulated class
    table: object                  # per-residue DataFrame

    @property
    def beats_baseline(self) -> bool:
        return bool(np.isfinite(self.agreement)
                    and self.agreement > self.baseline_agreement)

    def __repr__(self) -> str:
        verdict = "" if self.beats_baseline else "  [below majority baseline]"
        return (f"ShiftComparison(n={self.n_compared}, "
                f"MCC={self.mcc:.3f}, agreement={self.agreement:.1%} "
                f"vs {self.baseline_agreement:.1%} baseline, "
                f"balanced={self.balanced_accuracy:.1%}){verdict}")


def compare_to_dssp(csi_table, dssp_codes: np.ndarray, resids: np.ndarray,
                    offset: int = 0) -> ShiftComparison:
    """Match CSI labels against DSSP over the trajectory.

    ``dssp_codes`` is the per-frame DSSP array from ``mdtraj.compute_dssp`` in
    simplified mode; the modal assignment per residue is used. ``offset`` is
    added to structure residue numbers before matching BMRB numbering -- these
    two almost never agree, so check the overlap before trusting the result.
    """
    import pandas as pd

    modal = []
    for col in range(dssp_codes.shape[1]):
        vals, counts = np.unique(dssp_codes[:, col], return_counts=True)
        modal.append(vals[np.argmax(counts)])
    sim = pd.DataFrame({"resid": np.asarray(resids) + offset, "dssp": modal})
    # mdtraj simplified DSSP uses H / E / C already
    merged = csi_table.merge(sim, on="resid", how="inner")
    merged = merged[merged["csi"] != "-"]
    nan = float("nan")
    if merged.empty:
        return ShiftComparison(0, nan, nan, nan, nan, nan, nan, nan, None, merged)

    truth, pred = merged["csi"].values, merged["dssp"].values
    match = truth == pred

    def recall(label: str) -> float:
        sel = truth == label
        return float((pred[sel] == label).mean()) if sel.any() else nan

    # What a model containing no structural information would score: always
    # answer whichever class the NMR data labels most often.
    baseline = float(pd.Series(truth).value_counts(normalize=True).max())
    recalls = [recall(c) for c in "HEC"]
    confusion = pd.crosstab(pd.Series(truth, name="NMR"),
                            pd.Series(pred, name="simulated"))

    return ShiftComparison(
        n_compared=len(merged),
        agreement=float(match.mean()),
        baseline_agreement=baseline,
        balanced_accuracy=float(np.nanmean(recalls)),
        mcc=matthews_corrcoef(truth, pred),
        helix_agreement=recalls[0],
        strand_agreement=recalls[1],
        coil_agreement=recalls[2],
        confusion=confusion,
        table=merged,
    )


def threshold_scan(sec_shifts, dssp_codes: np.ndarray, resids: np.ndarray,
                   thresholds=(0.7, 1.0, 1.4, 1.75, 2.0, 2.5, 3.0, 3.5),
                   offset: int = 0):
    """How the comparison responds to the CSI threshold.

    The threshold is a judgement call, so report the sensitivity rather than a
    single number chosen by hand. If the conclusion only holds at one
    threshold, it is a property of that choice and not of the structure.
    """
    import pandas as pd

    rows = []
    for t in thresholds:
        cmp = compare_to_dssp(csi_assignment(sec_shifts, threshold=t),
                              dssp_codes, resids, offset=offset)
        rows.append({"threshold": t, "n": cmp.n_compared, "MCC": cmp.mcc,
                     "agreement": cmp.agreement,
                     "baseline": cmp.baseline_agreement,
                     "balanced": cmp.balanced_accuracy,
                     "beats_baseline": cmp.beats_baseline})
    return pd.DataFrame(rows)


def best_offset(csi_table, dssp_codes: np.ndarray, resids: np.ndarray,
                search: range = range(-20, 21)) -> tuple[int, ShiftComparison]:
    """Find the numbering offset that best aligns BMRB and structure residues.

    BMRB entries frequently number from the construct rather than the mature
    protein. Rather than assume, scan a window and take the best overlap.

    Selection is by **overlap count first**, with MCC only breaking ties.
    Numbering alignment is a question about which residues correspond, not
    about whether the model is right; scoring it by agreement would let an
    offset win because it happened to align well by chance, and would quietly
    bias the comparison it is supposed to set up.
    """
    best: tuple[int, ShiftComparison] | None = None
    best_key = None
    for off in search:
        cmp = compare_to_dssp(csi_table, dssp_codes, resids, offset=off)
        if cmp.n_compared == 0:
            continue
        key = (cmp.n_compared, cmp.mcc if np.isfinite(cmp.mcc) else -1.0)
        if best_key is None or key > best_key:
            best, best_key = (off, cmp), key
    if best is None:
        raise ValueError("no residue overlap at any offset in the search range")
    return best
