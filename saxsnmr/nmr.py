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

__all__ = ["secondary_shifts", "csi_assignment", "compare_to_dssp", "ShiftComparison"]

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

# A secondary shift smaller than this is not conformationally informative.
CSI_THRESHOLD = 0.7


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


def csi_assignment(sec_shifts):
    """Label each residue H (helix), E (strand) or C (coil) from its shifts."""
    import pandas as pd

    def label(row) -> str:
        v = row["d_diff"]
        if not np.isfinite(v):
            return "-"
        if v > CSI_THRESHOLD:
            return "H"
        if v < -CSI_THRESHOLD:
            return "E"
        return "C"

    out = sec_shifts.copy()
    out["csi"] = out.apply(label, axis=1)
    return out


@dataclass
class ShiftComparison:
    """Agreement between NMR-derived and simulated secondary structure."""

    n_compared: int
    agreement: float               # fraction of residues matching H/E/C
    helix_agreement: float
    strand_agreement: float
    table: object                  # per-residue DataFrame

    def __repr__(self) -> str:
        return (f"ShiftComparison(n={self.n_compared}, "
                f"agreement={self.agreement:.1%}, "
                f"helix={self.helix_agreement:.1%}, "
                f"strand={self.strand_agreement:.1%})")


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
    if merged.empty:
        return ShiftComparison(0, float("nan"), float("nan"), float("nan"), merged)

    match = merged["csi"] == merged["dssp"]

    def frac(label: str) -> float:
        sub = merged[merged["csi"] == label]
        return float((sub["csi"] == sub["dssp"]).mean()) if len(sub) else float("nan")

    return ShiftComparison(
        n_compared=int(len(merged)),
        agreement=float(match.mean()),
        helix_agreement=frac("H"),
        strand_agreement=frac("E"),
        table=merged,
    )


def best_offset(csi_table, dssp_codes: np.ndarray, resids: np.ndarray,
                search: range = range(-20, 21)) -> tuple[int, ShiftComparison]:
    """Find the numbering offset that best aligns BMRB and structure residues.

    BMRB entries frequently number from the construct rather than the mature
    protein. Rather than assume, scan a window and take the best overlap.
    """
    best: tuple[int, ShiftComparison] | None = None
    for off in search:
        cmp = compare_to_dssp(csi_table, dssp_codes, resids, offset=off)
        if cmp.n_compared == 0:
            continue
        score = cmp.n_compared * (cmp.agreement if np.isfinite(cmp.agreement) else 0)
        if best is None or score > best[1].n_compared * best[1].agreement:
            best = (off, cmp)
    if best is None:
        raise ValueError("no residue overlap at any offset in the search range")
    return best
