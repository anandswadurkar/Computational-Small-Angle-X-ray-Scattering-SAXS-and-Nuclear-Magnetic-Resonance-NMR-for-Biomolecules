"""Fetch and join experimental data for a biomolecule.

Sources, and what each is responsible for:

    SASBDB     experimental I(q), Rg, Dmax, sample/construct metadata
    PDBe/SIFTS UniProt -> ranked experimental structures (SASBDB rarely links one)
    RCSB       atomic coordinates
    AlphaFold  predicted structure, used only when no experimental one exists
    BMRB       assigned chemical shifts, mapped by UniProt

Two source quirks are handled here rather than left to callers:

1. SASBDB publishes q in either 1/A or 1/nm, declared in the entry's
   ``angular_unit`` field and *not* inferable from magnitude. Everything this
   module returns is normalised to 1/A, with the original unit recorded.
2. A UniProt accession identifies a protein, not a sample. Entries sharing an
   accession routinely differ in construct, tags, truncation and oligomeric
   state. ``constructs()`` surfaces those fields so a caller can refuse a join
   rather than average over incompatible samples.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import requests

__all__ = [
    "Fetcher",
    "SaxsCurve",
    "Pddf",
    "Construct",
    "StructureHit",
]

SASBDB = "https://www.sasbdb.org"
PDBE = "https://www.ebi.ac.uk/pdbe/api"
RCSB_FILES = "https://files.rcsb.org/download"
ALPHAFOLD = "https://alphafold.ebi.ac.uk/api/prediction"
BMRB_API = "https://api.bmrb.io/v2"

# SASBDB's declared angular units -> multiplier that converts q to 1/A.
_UNIT_TO_INV_ANGSTROM = {
    "1/a": 1.0,
    "1/angstrom": 1.0,
    "a-1": 1.0,
    "1/nm": 0.1,
    "nm-1": 0.1,
}


def _nm_to_A(value: float | None) -> float | None:
    """SASBDB reports Rg and Dmax in nm regardless of its angular unit."""
    return None if value is None else float(value) * 10.0


@dataclass
class Pddf:
    """An experimental pair-distance distribution, in Angstrom.

    This is GNOM's indirect Fourier transform of the measured intensity, not
    the raw data: its shape depends on the Dmax and regularisation the
    depositor chose. Treat it as a second opinion about the molecule's size,
    not as an independent measurement.
    """

    code: str
    r: np.ndarray
    p: np.ndarray
    error: np.ndarray
    source: str

    @property
    def dmax(self) -> float:
        """Where GNOM's P(r) was terminated -- a fitted choice, not a datum."""
        return float(self.r[-1])

    @property
    def rg(self) -> float:
        """Rg from the second moment, the same way `rg_from_pr` computes it."""
        total = self.p.sum()
        if total <= 0:
            return float("nan")
        return float(np.sqrt((self.r ** 2 * self.p).sum() / (2.0 * total)))

    def normalised(self) -> np.ndarray:
        peak = self.p.max()
        return self.p / peak if peak > 0 else self.p

    def __repr__(self) -> str:
        return (f"Pddf({self.code}, n={self.r.size}, "
                f"Dmax={self.dmax:.1f} A, Rg={self.rg:.2f} A)")


@dataclass
class SaxsCurve:
    """An experimental scattering curve, always in 1/A.

    The reported Rg and Dmax are converted to **Angstrom** to match q, which
    SASBDB does not do: it reports q in whichever unit `angular_unit` names but
    always gives Rg and Dmax in nanometres. Mixing the two silently rescales
    every length by ten, so the fields carry an explicit `_A` suffix. Anything
    still reading `curve.guinier_rg` will fail loudly rather than be wrong.
    """

    code: str
    q: np.ndarray
    intensity: np.ndarray
    sigma: np.ndarray | None
    original_unit: str
    guinier_rg_A: float | None = None
    guinier_rg_err_A: float | None = None
    pddf_rg_A: float | None = None
    pddf_dmax_A: float | None = None

    @property
    def q_range(self) -> tuple[float, float]:
        return float(self.q.min()), float(self.q.max())

    def trim(self, q_max: float = 0.5) -> "SaxsCurve":
        """Restrict to the range where an implicit hydration shell is valid.

        Beyond roughly 0.5 1/A the solvent's own structure contributes and a
        uniform shell stops being a good model. This project does not attempt
        that regime; see the README.
        """
        m = self.q <= q_max
        return SaxsCurve(
            code=self.code,
            q=self.q[m],
            intensity=self.intensity[m],
            sigma=None if self.sigma is None else self.sigma[m],
            original_unit=self.original_unit,
            guinier_rg_A=self.guinier_rg_A,
            guinier_rg_err_A=self.guinier_rg_err_A,
            pddf_rg_A=self.pddf_rg_A,
            pddf_dmax_A=self.pddf_dmax_A,
        )

    def __repr__(self) -> str:
        lo, hi = self.q_range
        return (f"SaxsCurve({self.code}, n={self.q.size}, "
                f"q={lo:.4f}-{hi:.3f} 1/A, from {self.original_unit})")


@dataclass
class Construct:
    """One molecule as actually deposited, with the fields a join must respect."""

    name: str | None
    molecular_type: str | None
    organism: str | None
    uniprot: str | None
    sequence: str | None

    @property
    def length(self) -> int | None:
        return len(self.sequence.replace("\n", "").strip()) if self.sequence else None


@dataclass
class StructureHit:
    pdb_id: str
    chain_id: str
    resolution: float | None
    method: str | None
    start: int | None
    end: int | None
    coverage: float | None = None


@dataclass
class Fetcher:
    """Cached access to the five sources. All network calls go through here."""

    cache_dir: Path = field(default_factory=lambda: Path(".cache"))
    timeout: int = 45

    def __post_init__(self) -> None:
        self.cache_dir = Path(self.cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._session = requests.Session()
        self._session.headers["User-Agent"] = "saxsnmr/0.1 (research use)"

    # ---------------------------------------------------------------- caching

    def _cached(self, key: str, url: str, binary: bool = False) -> Any:
        path = self.cache_dir / key
        if path.exists():
            return path.read_bytes() if binary else path.read_text(encoding="utf-8")
        r = self._session.get(url, timeout=self.timeout)
        r.raise_for_status()
        path.parent.mkdir(parents=True, exist_ok=True)
        if binary:
            path.write_bytes(r.content)
            return r.content
        # newline="" disables newline translation on write. Without it, text
        # that already contains CRLF -- which these servers send -- has its LF
        # translated again on Windows, producing CRLFCR... and a spurious blank
        # line between every real one when the file is read back. Parsers that
        # skip blanks survive that; anything that treats a blank line as the
        # end of a block silently truncates at the first row.
        path.write_text(r.text, encoding="utf-8", newline="")
        return r.text

    # ---------------------------------------------------------------- SASBDB

    def search_sasbdb(self, term: str) -> list[str]:
        """Find SASBDB entry codes for a free-text term.

        SASBDB exposes no working REST search endpoint -- every documented one
        returns HTTP 500 -- so this parses the HTML result page. It is the most
        fragile call in the module and will break if SASBDB restyles its site.
        """
        html = self._cached(f"sasbdb_search_{_slug(term)}.html",
                            f"{SASBDB}/search/?q={requests.utils.quote(term)}")
        return sorted(set(re.findall(r"SASD[A-Z]{2}\d", html)))

    def sasbdb_entry(self, code: str) -> dict:
        txt = self._cached(f"sasbdb_{code}.json",
                           f"{SASBDB}/rest-api/entry/summary/{code}/")
        return json.loads(txt)

    def constructs(self, code: str) -> list[Construct]:
        entry = self.sasbdb_entry(code)
        mols = (entry.get("experiment", {}).get("sample", {}) or {}).get("molecule", []) or []
        out = []
        for m in mols:
            out.append(Construct(
                name=m.get("name") or m.get("molecule_name"),
                molecular_type=m.get("molecular_type"),
                organism=m.get("organism"),
                uniprot=(m.get("uniprot_code") or None),
                sequence=m.get("sequence"),
            ))
        return out

    def sasbdb_curve(self, code: str) -> SaxsCurve:
        """Experimental curve, normalised to 1/A using the declared unit."""
        entry = self.sasbdb_entry(code)
        declared = (entry.get("angular_unit") or "").strip()
        scale = _UNIT_TO_INV_ANGSTROM.get(declared.lower())
        if scale is None:
            raise ValueError(
                f"{code}: unrecognised angular_unit {declared!r}. Refusing to guess "
                "-- inferring units from magnitude silently corrupts q by 10x."
            )

        text = self._cached(f"sasbdb_{code}.dat",
                            f"{SASBDB}/media/intensities_files/{code}.dat")
        q, i, s = [], [], []
        for line in text.splitlines():
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                vals = [float(p) for p in parts[:3]]
            except ValueError:
                continue
            if not np.isfinite(vals[0]) or vals[0] <= 0:
                continue
            q.append(vals[0])
            i.append(vals[1])
            s.append(vals[2] if len(vals) > 2 else np.nan)

        q = np.asarray(q) * scale
        sig = np.asarray(s)
        return SaxsCurve(
            code=code,
            q=q,
            intensity=np.asarray(i),
            sigma=None if np.all(np.isnan(sig)) else sig,
            original_unit=declared,
            # SASBDB reports these in nm whatever angular_unit says; see the
            # SaxsCurve docstring for why they are converted here.
            guinier_rg_A=_nm_to_A(entry.get("guinier_rg")),
            guinier_rg_err_A=_nm_to_A(entry.get("guinier_rg_error")),
            pddf_rg_A=_nm_to_A(entry.get("pddf_rg")),
            pddf_dmax_A=_nm_to_A(entry.get("pddf_dmax")),
        )

    def pddf(self, code: str) -> "Pddf":
        """Experimental P(r), parsed from the deposited GNOM output file.

        SASBDB links the raw `.out` file that ATSAS GNOM produced. Everything
        before the distance-distribution block is GNOM's own log -- settings,
        the regularised fit, perceptual criteria -- and is skipped. R is in
        Angstrom in these files and the final row sits exactly at Dmax.
        """
        entry = self.sasbdb_entry(code)
        url = entry.get("pddf_data")
        if not url:
            raise ValueError(f"{code}: no P(r) file deposited")

        text = self._cached(f"sasbdb_{code}_pofr.out", url)
        lines = text.splitlines()
        start = None
        for i, line in enumerate(lines):
            parts = line.split()
            if [p.upper() for p in parts[:3]] == ["R", "P(R)", "ERROR"]:
                start = i + 1
                break
        if start is None:
            raise ValueError(f"{code}: no 'R P(R) ERROR' block in the GNOM file")

        r, p, err = [], [], []
        for line in lines[start:]:
            parts = line.split()
            if not parts:
                continue           # blank lines separate nothing here
            if len(parts) != 3:
                if r:
                    break          # a real trailing section
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                if r:
                    break
                continue
            r.append(vals[0])
            p.append(vals[1])
            err.append(vals[2])

        if len(r) < 3:
            raise ValueError(f"{code}: parsed only {len(r)} P(r) points")
        return Pddf(code=code, r=np.asarray(r), p=np.asarray(p),
                    error=np.asarray(err), source=url)

    # ------------------------------------------------------------- structures

    def best_structures(self, uniprot: str) -> list[StructureHit]:
        """Experimental structures for an accession, best first (PDBe/SIFTS)."""
        try:
            txt = self._cached(f"pdbe_best_{uniprot}.json",
                               f"{PDBE}/mappings/best_structures/{uniprot}")
        except requests.HTTPError:
            return []
        hits = json.loads(txt).get(uniprot, [])
        return [StructureHit(
            pdb_id=h["pdb_id"], chain_id=h["chain_id"],
            resolution=h.get("resolution"),
            method=h.get("experimental_method"),
            start=h.get("unp_start") or h.get("start"),
            end=h.get("unp_end") or h.get("end"),
            coverage=h.get("coverage"),
        ) for h in hits]

    def pdb_file(self, pdb_id: str) -> Path:
        key = f"{pdb_id.lower()}.pdb"
        self._cached(key, f"{RCSB_FILES}/{pdb_id.lower()}.pdb", binary=True)
        return self.cache_dir / key

    def alphafold_file(self, uniprot: str) -> Path | None:
        try:
            txt = self._cached(f"af_{uniprot}.json", f"{ALPHAFOLD}/{uniprot}")
        except requests.HTTPError:
            return None
        meta = json.loads(txt)
        if not meta:
            return None
        url = meta[0]["pdbUrl"]
        key = f"AF_{uniprot}.pdb"
        self._cached(key, url, binary=True)
        return self.cache_dir / key

    def structure_for(self, uniprot: str) -> tuple[Path, str]:
        """Best available coordinates: experimental if any, else AlphaFold."""
        hits = self.best_structures(uniprot)
        if hits:
            top = hits[0]
            return self.pdb_file(top.pdb_id), f"PDB {top.pdb_id} chain {top.chain_id}"
        af = self.alphafold_file(uniprot)
        if af is None:
            raise LookupError(f"no structure available for {uniprot}")
        return af, f"AlphaFold {uniprot} (predicted)"

    # -------------------------------------------------------------------- NMR

    def bmrb_ids(self, uniprot: str) -> list[str]:
        txt = self._cached("bmrb_uniprot_map.json", f"{BMRB_API}/mappings/uniprot/bmrb")
        base = uniprot.split("-")[0].upper()
        for row in json.loads(txt):
            if row["uniprot_id"].split("-")[0].upper() == base:
                return list(row["bmrb_ids"])
        return []

    def bmrb_shifts(self, bmrb_id: str):
        """Assigned chemical shifts as a DataFrame, or None if the entry has none."""
        import pandas as pd

        txt = self._cached(
            f"bmrb_{bmrb_id}_cs.json",
            f"{BMRB_API}/entry/{bmrb_id}?saveframe_category=assigned_chemical_shifts")
        data = json.loads(txt).get(str(bmrb_id), {}).get("assigned_chemical_shifts", [])
        frames = []
        for sf in data:
            for loop in sf.get("loops", []):
                tags = [t.lstrip("_").split(".")[-1] for t in loop.get("tags", [])]
                rows = loop.get("data", [])
                if not rows or "Atom_ID" not in tags:
                    continue
                frames.append(pd.DataFrame(rows, columns=tags))
        if not frames:
            return None

        df = pd.concat(frames, ignore_index=True)
        keep = {"Comp_index_ID": "resid", "Seq_ID": "seq_id", "Comp_ID": "resname",
                "Atom_ID": "atom", "Atom_type": "element", "Val": "shift",
                "Val_err": "shift_err"}
        df = df[[c for c in keep if c in df.columns]].rename(columns=keep)
        for col in ("resid", "seq_id"):
            if col in df:
                df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
        for col in ("shift", "shift_err"):
            if col in df:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        return df.dropna(subset=["shift"])


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def summarise_constructs(fetcher: Fetcher, codes: Iterable[str]) -> "Any":
    """Side-by-side construct table for candidate entries.

    Intended to be read before joining anything. Differing sequence lengths or
    wildly differing Rg/Dmax across entries that share an accession mean the
    samples are not the same physical system.
    """
    import pandas as pd

    rows = []
    for code in codes:
        try:
            entry = fetcher.sasbdb_entry(code)
            for c in fetcher.constructs(code):
                rows.append({
                    "code": code,
                    "uniprot": c.uniprot,
                    "type": c.molecular_type,
                    "organism": c.organism,
                    "length": c.length,
                    # Angstrom, like every other length this package returns.
                    "guinier_rg_A": _nm_to_A(entry.get("guinier_rg")),
                    "pddf_dmax_A": _nm_to_A(entry.get("pddf_dmax")),
                    "angular_unit": entry.get("angular_unit"),
                    "type_of_curve": entry.get("type_of_curve"),
                })
        except Exception as exc:  # a dead entry should not kill the survey
            rows.append({"code": code, "uniprot": f"ERROR: {exc}"})
    return pd.DataFrame(rows)
