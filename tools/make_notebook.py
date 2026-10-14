"""Generate the validation notebooks, one per reference system.

The notebooks are build output; this file is the source. Editing cells in
Jupyter and saving works for a quick experiment, but the next run of this
script overwrites them, so anything worth keeping belongs here.

    python tools/make_notebook.py            # all systems
    python tools/make_notebook.py lysozyme   # just one

Every system runs the identical analysis. Only the molecule, the accessions
and the system-specific commentary differ, which is the point: a harness that
needed different handling per protein would not be a harness. Where a system
has its own trap -- a cofactor the force field cannot represent, a chain that
is not the first in the file -- it is named in that system's config rather
than worked around silently.

Hand-editing .ipynb JSON is error-prone and its diffs are unreadable, which is
the whole reason for generating it.
"""
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


# Shared by every notebook. Kept in one place because the two builders
# drifting apart is precisely what produced a kernel crash: the SANS
# notebook was written with its own import cell that put NumPy first.
BEFORE_YOU_RUN = """
## Before you run this

**1. Use a kernel from the `saxsnmr` environment.** On Windows, conda keeps its
compiled libraries in `<env>\\Library\\bin`, which only joins the search path
when the environment is activated. A Jupyter kernel launches `python.exe`
directly, so it normally misses them, and the first call into NumPy's linear
algebra then kills the kernel with no exception and no traceback — reported
only as *"the kernel crashed"*.

**2. `saxsnmr` must be imported before NumPy, SciPy, matplotlib, mdtraj or
pandas.** Importing it is what repairs PATH. Import NumPy first — or
matplotlib, which pulls NumPy in — and the compiled libraries bind against an
incomplete PATH before the repair runs, and the kernel dies on the first cell.

This ordering is load-bearing, not style. It is checked when these notebooks
are generated, and the first cell verifies the result. If you add imports to
the first cell, keep `import saxsnmr` above them.

**3. If you edit the package while a kernel is running, restart it.** Python
caches imported modules, so new names stay invisible and show up much later as
a confusing `ImportError`. To reload without losing your variables:

```python
import sys
for m in [k for k in sys.modules if k == "saxsnmr" or k.startswith("saxsnmr.")]:
    del sys.modules[m]
```

Purge the whole package, not just the top level — `importlib.reload` leaves the
submodules cached and mixes versions, which fails more confusingly than before.
"""


@dataclass
class System:
    """One reference molecule and everything that differs because of it."""

    index: str                  # notebook number prefix
    slug: str                   # file name stem after the number
    name: str                   # full name, for prose
    short: str                  # short name, used inline
    uniprot: str
    sasbdb_code: str
    search_term: str            # free-text SASBDB search
    q_max: float                # trim point, 1/A
    monomer_note: str           # how to spot a non-monomer in the survey
    structure_note: str         # what preparation has to cope with
    mutation: str               # PDBFixer notation
    mutation_label: str         # e.g. "E35Q"
    mutation_rationale: str     # why this residue
    why_this_system: str        # opening framing
    extra_caveats: list = field(default_factory=list)
    sans_code: str | None = None      # deposited SANS curve, if any
    sans_d2o: float = 0.0             # D2O fraction of that measurement
    n_replicas: int = 3
    production_ps: float = 5.0
    equilibration_ps: float = 2.0
    frame_interval_ps: float = 0.5

    @property
    def filename(self) -> str:
        return f"{self.index}_{self.slug}.ipynb"


SYSTEMS = [
    System(
        index="02", slug="phase1_lysozyme",
        name="hen egg white lysozyme", short="lysozyme",
        uniprot="P00698", sasbdb_code="SASDMJ2", search_term="lysozyme",
        q_max=0.50,
        monomer_note=(
            "Monomeric lysozyme is 129 residues with Rg ~14 A and Dmax ~42 A. "
            "Entries reporting a much larger Rg or Dmax are aggregates or "
            "fibrils rather than the folded monomer, and must not be averaged "
            "together with it. `SASDMJ4` is the clearest example here: it "
            "reports Rg ~303 A against the monomers' ~14 A."),
        structure_note=(
            "Lysozyme's best structure, `2vb1` at 0.65 A, needs all three: it "
            "resolves its own hydrogens, carries altloc codes A/B/C, and "
            "labels two different chains `A`."),
        mutation="GLU-35-GLN", mutation_label="E35Q",
        mutation_rationale=(
            "Glu35 is lysozyme's catalytic general acid. **E35Q** is a classic "
            "inactivating substitution: it removes the charge while keeping "
            "almost exactly the same shape. That makes it a good test of the "
            "method's limits, because a change that drastic for the enzyme's "
            "*function* should be nearly invisible to SAXS, which only sees "
            "the shape."),
        why_this_system=(
            "Lysozyme is the regression target because it is exhaustively "
            "characterised. If the harness cannot reproduce lysozyme, nothing "
            "else it reports is trustworthy."),
    ),
    System(
        index="03", slug="phase1_rnase_a",
        name="bovine pancreatic ribonuclease A", short="RNase A",
        uniprot="P61823", sasbdb_code="SASDPF7", search_term="ribonuclease",
        q_max=0.50,
        monomer_note=(
            "Monomeric RNase A is 124 residues with Rg ~15 A and Dmax ~46 A. "
            "The entries here were collected by SEC-SAXS, where the sample is "
            "separated immediately before measurement, so aggregates are far "
            "less of a concern than in a batch experiment -- the four entries "
            "agree on Rg to within 0.2 A."),
        structure_note=(
            "RNase A's best structure, `6etk` at 0.85 A, is comparable to "
            "lysozyme's reference: high enough resolution to carry its own "
            "hydrogens and alternate conformations. It has four disulfides, "
            "which the force field handles but which make the fold stiff."),
        mutation="HIS-12-ALA", mutation_label="H12A",
        mutation_rationale=(
            "His12 and His119 are RNase A's catalytic pair, acting as general "
            "base and general acid. **H12A** abolishes activity. Like E35Q in "
            "lysozyme it is a chemical change rather than a shape change, so "
            "SAXS should not see it -- but unlike E35Q it removes a side chain "
            "rather than swapping one, so PDBFixer's rebuilt rotamer matters "
            "even less."),
        why_this_system=(
            "RNase A is the closest available analogue to lysozyme: a small, "
            "rigid, exhaustively studied enzyme with a sub-Angstrom structure. "
            "It exists here to answer whether the lysozyme result is a "
            "property of lysozyme or of the method. A conclusion that holds "
            "for one protein and not the other is a different kind of finding "
            "from one that holds for both."),
        extra_caveats=[
            "The deposited curve is in 1/nm, not 1/A. The loader converts it "
            "from the `angular_unit` field, and the deposited P(r) is "
            "converted the same way -- GNOM writes r in whatever length unit "
            "matches the input q, which is a silent factor of ten if assumed.",
        ],
    ),
    System(
        index="04", slug="phase1_cytochrome_c",
        name="equine heart cytochrome c", short="cytochrome c",
        uniprot="P00004", sasbdb_code="SASDAB2", search_term="cytochrome c",
        q_max=0.50,
        monomer_note=(
            "Monomeric cytochrome c is ~104 residues with Rg ~13 A and "
            "Dmax ~37 A. It is the smallest of the three reference systems, "
            "so its simulations are correspondingly quicker."),
        structure_note=(
            "Cytochrome c's structure is solved as part of a larger assembly, "
            "so the chain of interest is **not the first in the file**. The "
            "chain identifier comes from the PDBe/SIFTS mapping rather than "
            "being assumed, and `sanitise_pdb` has to discard the TER and "
            "SSBOND records belonging to the chains it drops."),
        mutation="MET-80-ALA", mutation_label="M80A",
        mutation_rationale=(
            "Met80 provides the sixth coordination bond to the heme iron, and "
            "**M80A** is the classic substitution for opening that site. Note "
            "the asymmetry with the other two systems: here the mutation's "
            "whole significance is to a cofactor this simulation does not "
            "contain (see the caveat below), so it is being used purely as a "
            "shape perturbation."),
        why_this_system=(
            "Cytochrome c is a different fold class from the other two and it "
            "carries a cofactor, which makes it the most honest of the three "
            "about what this pipeline cannot do. It is included because a "
            "harness that only ever sees cofactor-free globular enzymes will "
            "not warn you when that assumption breaks."),
        extra_caveats=[
            "**The heme is not simulated.** amber14 has no parameters for "
            "Fe-protoporphyrin IX, and `prepare` strips heterogens, so this "
            "is apo-cytochrome c. The crystal structure still fits the "
            "experimental curve well because the heme is a small fraction of "
            "the total scattering, but the simulated protein is missing the "
            "group its fold is built around. Treat any dynamics here as "
            "indicative only -- over longer timescales apo-cytochrome c is "
            "known to be a molten globule, and a longer run would show that "
            "rather than the native state the experiment measured.",
        ],
    ),
]


class Notebook:
    """Accumulates cells for one notebook."""

    def __init__(self, system: System):
        self.s = system
        self.cells = []
        self._n = 0

    def _id(self):
        self._n += 1
        return f"cell{self._n:03d}"

    def _fill(self, text: str) -> str:
        s = self.s
        return (text
                .replace("__UNIPROT__", s.uniprot)
                .replace("__CODE__", s.sasbdb_code)
                .replace("__TERM__", s.search_term)
                .replace("__NAME__", s.name)
                .replace("__SHORT__", s.short)
                .replace("__QMAX__", f"{s.q_max:.2f}")
                .replace("__MUT__", s.mutation)
                .replace("__MUTLABEL__", s.mutation_label)
                .replace("__NREP__", str(s.n_replicas))
                .replace("__PRODPS__", f"{s.production_ps}")
                .replace("__EQPS__", f"{s.equilibration_ps}")
                .replace("__FRAMEPS__", f"{s.frame_interval_ps}"))

    def md(self, text):
        self.cells.append({
            "cell_type": "markdown", "id": self._id(), "metadata": {},
            "source": self._fill(text).strip("\n").splitlines(keepends=True)})

    def code(self, text):
        self.cells.append({
            "cell_type": "code", "id": self._id(), "execution_count": None,
            "metadata": {}, "outputs": [],
            "source": self._fill(text).strip("\n").splitlines(keepends=True)})

    def write(self):
        # saxsnmr must be imported before numpy or matplotlib, in the notebook
        # as a whole. A Jupyter kernel runs without the conda environment
        # activated, so importing saxsnmr is what puts <env>\Library\bin on
        # PATH; bind BLAS before that and the process aborts with no traceback
        # on the first cell. This shipped once, and headless testing cannot
        # catch it because `conda run` activates the environment and makes the
        # repair a no-op.
        # Order is (cell index, line index within that cell). Comparing the
        # line text instead would order two imports in the same cell
        # alphabetically, which is meaningless and wrong -- the first version
        # of this check did exactly that and reported a false positive.
        first_sci, first_saxsnmr = None, None
        for n, cell in enumerate(self.cells, 1):
            if cell["cell_type"] != "code":
                continue
            for ln, line in enumerate(cell["source"]):
                stripped = line.strip()
                if not stripped.startswith(("import ", "from ")):
                    continue
                where = (n, ln)
                if first_saxsnmr is None and "saxsnmr" in stripped:
                    first_saxsnmr = (where, stripped)
                if first_sci is None and any(
                        stripped.startswith(p) for p in
                        ("import numpy", "import matplotlib", "from numpy",
                         "from matplotlib", "import scipy", "from scipy",
                         "import mdtraj", "from mdtraj", "import pandas")):
                    first_sci = (where, stripped)
        if first_sci and (first_saxsnmr is None
                          or first_sci[0] < first_saxsnmr[0]):
            raise RuntimeError(
                f"{self.s.filename}: {first_sci[1]!r} (cell {first_sci[0][0]}, "
                f"line {first_sci[0][1] + 1}) comes before importing saxsnmr"
                + (f" ({first_saxsnmr[1]!r}, cell {first_saxsnmr[0][0]}, "
                   f"line {first_saxsnmr[0][1] + 1})"
                   if first_saxsnmr else " (never imported)")
                + ".\nsaxsnmr repairs PATH for the conda environment's "
                  "compiled libraries. Importing NumPy first binds BLAS "
                  "before the repair and the kernel dies with no traceback.")

        # Compile every code cell before writing. The cell sources are built
        # by string substitution into templates, and a mangled escape yields
        # a notebook that looks fine until someone runs it -- which has
        # happened. Failing here costs a second; failing in a notebook costs
        # however long it took to reach that cell.
        for n, cell in enumerate(self.cells, 1):
            if cell["cell_type"] != "code":
                continue
            src = "".join(cell["source"])
            try:
                compile(src, f"{self.s.filename} cell {n}", "exec")
            except SyntaxError as exc:
                raise SyntaxError(
                    f"{self.s.filename} cell {n} does not compile: {exc}\n"
                    f"--- cell source ---\n{src}") from exc

        nb = {"cells": self.cells,
              "metadata": {
                  "kernelspec": {"display_name": "Python (saxsnmr)",
                                 "language": "python", "name": "saxsnmr"},
                  "language_info": {"name": "python", "version": "3.11"}},
              "nbformat": 4, "nbformat_minor": 5}
        out = REPO / "notebooks" / self.s.filename
        out.write_text(json.dumps(nb, indent=1), encoding="utf-8")
        return out, len(self.cells)


def build(s: System) -> Notebook:
    nb = Notebook(s)
    md, code = nb.md, nb.code
    runtime_min = int(6 + s.n_replicas * 2 * 3.5)

    md(f"""
# Validating short implicit-solvent MD against SAXS and NMR — {s.name}

Reference system: **{s.name}**, UniProt `{s.uniprot}`.

This notebook runs the full loop end to end:

1. survey the candidate SASBDB entries and **check they describe the same molecule**
2. load one experimental scattering curve, normalised to 1/A
3. **check the curve against itself** — recompute Guinier Rg from the raw data,
   and look at the Kratky plot before assuming the particle is folded
4. source atomic coordinates through UniProt
5. prepare and simulate in implicit solvent
6. compute an ensemble SAXS profile and fit it
7. compare the distance distribution P(r) with the deposited GNOM transform
8. compare NMR chemical shifts to the simulated secondary structure, against
   two controls
9. run a point mutant **with replicates**, to see whether the method can
   resolve the change at all

{s.why_this_system}

**Expect roughly {runtime_min}–{runtime_min + 10} minutes end to end on a laptop CPU**,
almost all of it in section 7, which runs {s.n_replicas * 2} short simulations to get
replicates. Sections 1–6 take a few minutes. Set `N_REPLICAS = 1` there if you
only want to see the pipeline work — but read what that section says about why
one run per group is not a control before trusting any number it then prints.
""")

    md(BEFORE_YOU_RUN)

    code("""
# Environment check. Run this first.
#
# Importing saxsnmr before anything else matters: on Windows it repairs PATH
# so the conda environment's compiled libraries can be found. Without that,
# the first call into NumPy's linear algebra kills the kernel outright -- no
# exception, no traceback, just "the kernel crashed".
import sys
from pathlib import Path

sys.path.insert(0, str(Path.cwd().parent))
import saxsnmr

# A kernel that imported saxsnmr before the package was last edited keeps the
# old module object, and new names are simply absent. That surfaces much
# later as an opaque ImportError in whichever cell first needs one, so check
# the API surface here instead.
REQUIRED = ["Fetcher", "prepare", "run_md", "ensemble_fit", "guinier",
            "pair_distribution", "threshold_scan", "shift_dispersion",
            "sans_tables", "sans_profile", "match_point_d2o"]
missing = [n for n in REQUIRED if not hasattr(saxsnmr, n)]
if missing:
    raise ImportError(
        f"the imported saxsnmr is out of date -- missing {missing}.\\n"
        "The kernel cached it before the package was edited. Either restart\\n"
        "the kernel, or run this to reload without losing your variables:\\n\\n"
        "    import sys\\n"
        "    for m in [k for k in sys.modules\\n"
        "              if k == 'saxsnmr' or k.startswith('saxsnmr.')]:\\n"
        "        del sys.modules[m]\\n\\n"
        "then re-run this cell. Reloading only the top-level package is not\\n"
        "enough: the submodules stay cached and you get a mix of versions.")

print("python:", sys.executable)
print("saxsnmr API complete:", len(REQUIRED), "entry points present")
if saxsnmr.DLL_PATH_ADDED:
    print("repaired PATH with", len(saxsnmr.DLL_PATH_ADDED), "library directories")
else:
    print("PATH already complete (environment was activated)")

checks = {
    "numpy linalg":   "import numpy as np; np.linalg.lstsq(np.eye(3), np.ones(3), rcond=None)",
    "scipy optimize": "from scipy.optimize import least_squares; least_squares(lambda p: p - 3.0, [0.0])",
    # nnls is a separate compiled entry point used by the SANS fit. Checking
    # least_squares does not exercise it, and anything that aborts the
    # process here would otherwise do so mid-notebook with no traceback.
    "scipy nnls":     ("import numpy as np; from scipy.optimize import nnls; "
                       "nnls(np.eye(3), np.ones(3))"),
    "mdtraj":         "import mdtraj",
    "openmm import":  "import openmm",
    "openmm runtime": ("import openmm as mm; s = mm.System(); s.addParticle(1.0); "
                       "mm.Context(s, mm.VerletIntegrator(0.001), "
                       "mm.Platform.getPlatformByName('CPU'))"),
    "pdbfixer":       "import pdbfixer",
}

failed = []
for name, snippet in checks.items():
    try:
        exec(snippet, {})
        print(f"  OK    {name}")
    except Exception as exc:
        print(f"  FAIL  {name}: {type(exc).__name__}: {exc}")
        failed.append(name)

if failed:
    print()
    print("Unusable:", ", ".join(failed))
    print("Select the 'Python (saxsnmr)' kernel, listed under Jupyter Kernels.")
else:
    import openmm
    print("\\nenvironment OK -- openmm", openmm.__version__)
""")

    code("""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path.cwd().parent))

from saxsnmr import (Fetcher, prepare, run_md, benchmark, ensemble_fit,
                     frame_tables, radius_of_gyration, summarise_constructs,
                     load_trajectory, guinier, pair_distribution,
                     rg_from_pr, dmax_from_pr)
from saxsnmr.saxs import ScatteringTables, fit
from saxsnmr.nmr import (secondary_shifts, csi_assignment, best_offset,
                         threshold_scan)

# Simulation output is per system; the download cache is shared, so running a
# second notebook does not re-fetch what the first already has.
WORK = Path("work") / "__SLUG__"
WORK.mkdir(parents=True, exist_ok=True)
fetcher = Fetcher(cache_dir=Path("work") / "cache")
UNIPROT = "__UNIPROT__"
""".replace("__SLUG__", s.slug.replace("phase1_", "")))

    md("""
## 1. Which entries actually describe the same molecule?

A shared UniProt accession identifies a *protein*, not a *sample*. Entries
differ in construct, tags, truncation and oligomeric state.

Read this table before trusting any join. Look for disagreement in
`length`, and for `guinier_rg_A` / `pddf_dmax_A` that do not match the rest —
those entries are a different physical system. Note also that `angular_unit`
varies between entries; the loader normalises it, but never assume.

**On units.** SASBDB reports q in whichever unit `angular_unit` names, but
always reports Rg and Dmax in nanometres regardless. Mixing the two rescales
every length by ten. This package converts both to Angstrom on the way in and
marks the fields with an `_A` suffix, so the mistake cannot be made silently.
""")

    code("""
codes = fetcher.search_sasbdb("__TERM__")
print(f"{len(codes)} candidate entries")

survey = summarise_constructs(fetcher, codes[:12])
survey[["code", "uniprot", "length", "guinier_rg_A", "pddf_dmax_A",
        "angular_unit", "type_of_curve"]]
""")

    md(f"""
### Reading the survey

{s.monomer_note}
""")

    code("""
CODE = "__CODE__"

curve = fetcher.sasbdb_curve(CODE)
print(curve)
print(f"SASBDB Guinier Rg = {curve.guinier_rg_A:.2f} +/- {curve.guinier_rg_err_A or float('nan'):.2f} A")
print(f"SASBDB P(r):  Rg = {curve.pddf_rg_A:.2f} A, Dmax = {curve.pddf_dmax_A:.1f} A")

# Restrict to where an implicit hydration shell is valid (see README).
curve = curve.trim(q_max=__QMAX__)
print("after trim:", curve)

fig, ax = plt.subplots(figsize=(6, 4))
ax.errorbar(curve.q, curve.intensity, yerr=curve.sigma, fmt=".", ms=3,
            lw=0.5, alpha=0.7)
ax.set_yscale("log")
ax.set_xlabel("q  (1/A)")
ax.set_ylabel("I(q)  (arb.)")
ax.set_title(f"{CODE} experimental")
plt.tight_layout()
""")

    md("""
## 1b. Is this curve what it claims to be?

Before using an experimental curve to judge a model, check the curve. Two
standard diagnostics do most of the work, and neither depends on having a
structure.

**Guinier.** At low q, `ln I(q) = ln I0 - Rg^2 q^2 / 3`. The approximation
holds only while `q*Rg` is small — conventionally below about 1.3 for a
globular particle. That criterion is circular, since the valid range depends
on the Rg being measured, so `guinier()` extends the window point by point and
stops when the fitted Rg puts the last point past the limit. A curve with no
valid window is flagged rather than fitted: aggregation shows up here as an
upturn at the lowest q, and it would otherwise be reported as a large Rg with
no indication anything is wrong.

This recomputes Rg from the raw data rather than trusting the deposited value.
Agreement between the two is a check on the whole loading path — units
included.
""")

    code("""
g = guinier(curve.q, curve.intensity, curve.sigma)
print(g)
print(f"  fitted {g.n_points} points, q = {g.q[0]:.4f} - {g.q[-1]:.4f} 1/A")
print(f"  SASBDB reports {curve.guinier_rg_A:.2f} A")
print(f"  difference: {abs(g.rg - curve.guinier_rg_A):.2f} A")
if g.note:
    print("  note:", g.note)
""")

    md("""
Note the error bars. This fit's formal error is far smaller than the one
SASBDB quotes, because it is the covariance of a straight line that assumes
the deposited sigmas are exactly right. The real uncertainty is dominated by
the choice of fitting range, which no covariance captures. `reduced_chi2` is
the more useful number: near 1 means the data really is a straight line over
the window; much above 1 means curvature, and therefore that the window is too
wide or the sample is not monodisperse.
""")

    code("""
fig, (axg, axk) = plt.subplots(1, 2, figsize=(11, 4))

# Guinier: a straight line here, curvature means trouble.
axg.errorbar(g.q ** 2, g.ln_intensity, fmt=".", ms=4, label="experiment")
axg.plot(g.q ** 2, g.ln_model, lw=1.5, color="crimson",
         label=f"fit: Rg = {g.rg:.2f} A")
axg.set_xlabel("q^2  (1/A^2)")
axg.set_ylabel("ln I(q)")
axg.set_title(f"Guinier region (qRg {g.q_rg_min:.2f}-{g.q_rg_max:.2f})")
axg.legend()

# Kratky: a folded globular particle gives a bell with a clear maximum that
# returns towards the axis. A curve that rises and keeps rising indicates a
# disordered or unfolded chain -- worth knowing before blaming the force field.
axk.plot(curve.q, curve.q ** 2 * curve.intensity, lw=1)
axk.axvline(1.73 / g.rg, color="grey", ls=":",
            label="qRg = 1.73 (expected peak)")
axk.set_xlabel("q  (1/A)")
axk.set_ylabel("q^2 I(q)")
axk.set_title("Kratky -- folded or not?")
axk.legend(fontsize=8)
plt.tight_layout()
""")

    md("""
For a compact, well-folded particle the Kratky peak sits near `qRg = 1.73`.
If a target does not show that, the premise that one short trajectory
represents its solution ensemble is already in doubt, and a poor chi-squared
later would be a symptom rather than a cause.

## 2. Coordinates

SASBDB links a structure for only about 15% of entries, so coordinates come
from PDBe/SIFTS keyed on UniProt, falling back to AlphaFold.

Preparation strips pre-existing hydrogens, resolves alternate conformations and
selects a single chain. """ + s.structure_note + """

The chain identifier is taken from the SIFTS mapping rather than assumed to be
`A`. For a protein solved inside a larger assembly it very often is not.
""")

    code("""
hits = fetcher.best_structures(UNIPROT)[:5]
for h in hits:
    print(f"  {h.pdb_id} chain {h.chain_id}  {h.resolution} A  {h.method}")

CHAIN = hits[0].chain_id
structure, provenance = fetcher.structure_for(UNIPROT)
print("\\nusing:", provenance, "chain", CHAIN)

prep = prepare(structure, WORK / "prepared.pdb", chain=CHAIN)
print(f"{prep.n_residues} residues, {prep.n_atoms} atoms")
print(prep.sequence)
""")

    md("""
## 3. How fast can this machine actually simulate?

Run this once. Throughput varies by more than an order of magnitude between
platforms, and every decision about trajectory length depends on the real number
rather than an estimate.

**Do not use OpenCL on integrated Intel graphics.** `python -m
openmm.testInstallation` on the development machine reported forces differing
from the CPU platform by 0.13, where CPU and Reference agree to 6e-06. A
platform that computes the wrong forces will still happily produce a
trajectory. Check your own output before trusting any platform other than CPU.
""")

    code("""
rates = benchmark(prep, ps=2.0, platforms=["CPU"])   # OpenCL is wrong here; Reference is very slow
for name, ns_day in rates.items():
    print(f"  {name:<10} {ns_day:8.2f} ns/day" if np.isfinite(ns_day)
          else f"  {name:<10}      unavailable")
""")

    md("""
## 4. Short implicit-solvent MD

Implicit solvent is what puts this on a laptop. It also carries a known bias,
and Rg is exactly what SAXS measures most sensitively — so treat the fit below
as a check that the workflow is sound, not as a verdict on the force field.

**One setting worth knowing about.** `run_md` leaves Debye screening at
OpenMM's default, which is zero: the simulation has unscreened intramolecular
electrostatics regardless of the salt the experiment was measured in. Testing
this on lysozyme showed it does *not* explain the expansion seen there — but
it is a mismatch with the experiment either way, and the buffer conditions are
in the SASBDB entry if you want to check them for this system.
""")

    code("""
# run_md's tuned defaults give roughly 4-8 ns/day for a system this size on a
# laptop CPU -- several times what the obvious settings give, mostly from a
# shorter nonbonded cutoff. Raise production_ps once the benchmark above tells
# you your rate.
result = run_md(prep, WORK, production_ps=__PRODPS__, equilibration_ps=__EQPS__,
                frame_interval_ps=__FRAMEPS__, platform="CPU")
print(result)
""")

    code("""
import mdtraj as mdt

# load_trajectory pairs the DCD with its topology and converts units.
traj = load_trajectory(result.trajectory, top=result.topology)
rg = radius_of_gyration(traj)

fig, ax = plt.subplots(figsize=(6, 3))
ax.plot(np.arange(len(rg)) * __FRAMEPS__, rg, lw=1)
ax.axhline(g.rg, color="crimson", ls="--",
           label=f"experimental Guinier Rg = {g.rg:.2f} A")
ax.legend()
ax.set_xlabel("time (ps)")
ax.set_ylabel("Rg (A)")
ax.set_title("Compaction check")
plt.tight_layout()
print(f"simulated Rg = {rg.mean():.2f} +/- {rg.std():.2f} A")
""")

    md("""
### 4b. Is this trajectory long enough to quote a mean?

Successive MD frames are not independent. Treating them as if they were makes
the standard error too small by roughly the square root of the correlation
time, which is how a run that has not converged reports a tight, confident and
wrong number.

`block_average` corrects for that, and `assess_convergence` checks two
distinct failures: a series can **drift**, in which case its mean describes no
particular state, or it can be stationary but so **short** that the mean is
unresolved. The fixes differ — longer equilibration against a longer
production run — so they are reported separately.

A 5 ps run sampled every 0.5 ps gives ten frames. Expect that to fail.
""")

    code("""
from saxsnmr import assess_convergence, block_average

stats = block_average(rg)
print(stats)
print(f"  naive standard error   : {stats.naive_sem:.4f} A")
print(f"  correlation-corrected  : {stats.sem:.4f} A  "
      f"({stats.inflation:.1f}x larger)")
print(f"  effective samples      : {stats.n_effective:.1f} of {stats.n_samples}")

report = assess_convergence(rg, "Rg")
print()
print(report)
# The two effective-sample counts differ on purpose: block_average above used
# every frame, while assess_convergence first discards what looks like
# equilibration and reports on the remainder.
print(f"  discarded as equilibration: {report.equilibration_index} of "
      f"{len(rg)} frames")
print(f"  effective samples in what remains: "
      f"{report.production.n_effective:.1f}")
if not report.converged:
    print("\\nThis is the expected answer at this trajectory length. The")
    print("systematic expansion is still visible -- it is several times the")
    print("scatter -- but the mean Rg should be quoted with the corrected")
    print("error, not the naive one.")
""")

    md("""
## 5. SAXS: ensemble profile and fit

The experiment measures one ensemble-averaged curve, so the model must be
averaged over frames *before* fitting rather than fitting frames individually.

`c1` (excluded volume) and `c2` (hydration contrast) are free parameters. They
are reported alongside chi-squared deliberately — with two parameters floating,
a good fit alone does not establish much.
""")

    code("""
ens = ensemble_fit(result, curve, stride=1)
print(ens)

static = load_trajectory(prep.path)
st = frame_tables(static, 0)
single = fit(st, curve.q, curve.intensity, curve.sigma)
print("single crystal structure:", single)
""")

    code("""
fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.5, 6), sharex=True,
                              gridspec_kw={"height_ratios": [3, 1]})
ax.errorbar(curve.q, curve.intensity, yerr=ens.sigma, fmt=".", ms=3, lw=0.5,
            alpha=0.5, label="experiment")
ax.plot(curve.q, single.model, lw=1.2, label=f"crystal structure, chi2={single.chi2:.2f}")
ax.plot(curve.q, ens.model, lw=1.5, label=f"MD ensemble, chi2={ens.chi2:.2f}")
ax.set_yscale("log")
ax.set_ylabel("I(q)")
ax.legend()
ax.set_title(f"{CODE} vs model")

axr.axhline(0, color="k", lw=0.5)
axr.plot(curve.q, (ens.model - curve.intensity) / ens.sigma, lw=1)
axr.set_xlabel("q  (1/A)")
axr.set_ylabel("residual / sigma")
plt.tight_layout()
""")

    md("""
## 5b. The same comparison in real space

Chi-squared is computed in reciprocal space, where a scale factor and two free
parameters are floating. P(r) — the distribution of interatomic distances —
gives an independent view with none of those adjustments. A model can be made
to fit I(q) reasonably while having the wrong shape; it is harder to hide that
in P(r).

Two caveats keep this honest. The experimental P(r) is not data: it is GNOM's
indirect Fourier transform of the intensity, and its shape depends on the Dmax
and regularisation the depositor chose. And the model's P(r) is weighted by
electron count with no hydration shell, where GNOM's reflects contrast against
solvent. Compare the *shape and extent*, not the fine detail.
""")

    code("""
exp_pr = fetcher.pddf(CODE)
print("experiment:", exp_pr)

# Every fit carries the geometry it was computed from, so the ensemble-averaged
# tables are available here without rebuilding them.
r_x, p_x = pair_distribution(single.tables)
r_m, p_m = pair_distribution(ens.tables)

print(f"\\n{'':22s} {'Rg (A)':>8s}  {'Dmax (A)':>9s}")
print(f"{'experiment (GNOM)':22s} {exp_pr.rg:8.2f}  {exp_pr.dmax:9.1f}")
print(f"{'crystal structure':22s} {rg_from_pr(r_x, p_x):8.2f}  {dmax_from_pr(r_x, p_x):9.1f}")
print(f"{'MD ensemble':22s} {rg_from_pr(r_m, p_m):8.2f}  {dmax_from_pr(r_m, p_m):9.1f}")
""")

    code("""
fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(exp_pr.r, exp_pr.normalised(), lw=2, color="k", label="experiment (GNOM)")
ax.fill_between(exp_pr.r,
                exp_pr.normalised() - exp_pr.error / exp_pr.p.max(),
                exp_pr.normalised() + exp_pr.error / exp_pr.p.max(),
                color="k", alpha=0.15, lw=0)
ax.plot(r_x, p_x, lw=1.3, label="crystal structure")
ax.plot(r_m, p_m, lw=1.3, label="MD ensemble")
ax.axvline(exp_pr.dmax, color="grey", ls=":", lw=1,
           label=f"experimental Dmax = {exp_pr.dmax:.0f} A")
ax.set_xlim(0, exp_pr.dmax * 1.25)
ax.set_xlabel("r  (A)")
ax.set_ylabel("P(r), normalised")
ax.set_title("Distance distribution")
ax.legend(fontsize=8)
plt.tight_layout()
""")

    md("""
A model whose P(r) extends past the experimental Dmax is too large — and
because SAXS is most sensitive to exactly that, it is the first thing to check
when chi-squared is poor.

## 5c. Neutron scattering — a prediction, not a validation
""" + (f"""
**There is no deposited SANS curve for {s.short}.** SASBDB has none for this
molecule, so everything below is a forward prediction with nothing to check it
against. It is included because the contrast behaviour is informative in its
own right, and because the same code *is* validated against experiment in
`05_sans_validation.ipynb` — read that one before trusting these numbers.
""" if not s.sans_code else "") + """
Neutrons scatter from nuclei, so the scattering length carries no q
dependence — the calculation is simpler than the X-ray case, and reuses the
same distance histograms. What makes it physically different is hydrogen:

- **Hydrogen scatters negatively** (b = −3.74 fm) where deuterium is strongly
  positive (+6.67 fm). Nothing in X-ray scattering behaves like this.
- **Labile hydrogens exchange with solvent.** Those on N, O and S swap with the
  buffer and become deuterium in D₂O; carbon-bound ones do not. The two
  classes are therefore tracked separately, from the bond graph rather than
  from atom names.
- **Hydrogen scatters incoherently** — 80.3 barn against 1.76 coherent — which
  contributes a large flat background carrying no structural information.
""")

    code("""
from saxsnmr import (sans_tables, sans_profile, sans_fit, sans_elements,
                    solvent_sld, match_point_d2o)
from saxsnmr.sans import SLD_H2O, SLD_D2O

labels = sans_elements(static.topology)
n_h = labels.count("H")
n_hx = labels.count("Hx")
print(f"hydrogens: {n_h} carbon-bound (fixed), {n_hx} on N/O/S (exchangeable)")
print(f"solvent SLD: H2O {SLD_H2O:+.4f}, D2O {SLD_D2O:+.4f} fm/A^3")

sans_tab = sans_tables(static, 0)
match = match_point_d2o(sans_tab)
print(f"\\ncontrast match point for this molecule: {match:.1%} D2O")
print("Literature quotes 40-45% for a typical protein. This runs slightly")
print("higher because complete exchange of labile positions is assumed, which")
print("overstates deuteration -- buried amides exchange slowly or not at all.")
""")

    md("""
At the match point the molecule's mean scattering length density equals the
solvent's, the coherent signal very nearly vanishes, and what remains is
dominated by the hydration shell — which this implementation does not model.
**That is where these numbers are least trustworthy**, not merely noisiest.
""")

    code("""
q_n = np.linspace(0.01, 0.30, 80)
fig, (axc, axi) = plt.subplots(1, 2, figsize=(11, 4))

for x in (0.0, 0.25, 0.42, 0.70, 1.0):
    axc.plot(q_n, sans_profile(sans_tab, q_n, x), lw=1.2, label=f"{x:.0%} D2O")
axc.axhline(0, color="k", lw=0.5)
axc.set_yscale("symlog")
axc.set_xlabel("q  (1/A)")
axc.set_ylabel("I(q), arbitrary units")
axc.set_title("Predicted contrast series")
axc.legend(fontsize=8)

frac = np.linspace(0, 1, 60)
i0 = [sans_profile(sans_tab, np.array([0.01]), x)[0] for x in frac]
axi.plot(frac, np.sqrt(np.abs(i0)), lw=1.5)
axi.axvline(match, color="crimson", ls="--", label=f"match point {match:.0%}")
axi.set_xlabel("D2O fraction")
axi.set_ylabel("sqrt(|I(0)|)  ~ |contrast|")
axi.set_title("Contrast vanishes at the match point")
axi.legend(fontsize=8)
plt.tight_layout()
""")

    md("""
The right panel is the signature of a contrast-variation experiment: the
square root of I(0) falls linearly to zero at the match point and rises again
beyond it. Reproducing that V shape is the weakest possible check on the
calculator — it only confirms the contrast arithmetic, not the shape of the
molecule.

## 6. NMR: chemical shifts against simulated secondary structure

Predicting raw shifts needs SPARTA+, SHIFTX2 or UCBShift, none of which builds
on Windows. The portable route is secondary chemical shifts: subtracting
random-coil values leaves a signal whose sign reports helix versus strand, which
can be compared with DSSP over the trajectory.

**This section needs controls more than any other.** An agreement percentage
on its own is uninterpretable: secondary structure classes are nowhere near
equally populated, so a model that always answers "helix" can score well while
containing no information at all. Two reference points are computed alongside
the simulation — the crystal structure, which is as close to correct as a
model gets, and a majority-class predictor, which is as close to useless as a
model gets. The MD result only means something between them.
""")

    md("""
### Choosing a BMRB entry

Several entries usually exist for one protein, and they are not
interchangeable. Taking the first with enough CA shifts is not good enough:
for RNase A that picks an entry whose secondary shifts have a standard
deviation of 0.59 ppm where five other entries give 3.6–5.0 ppm. An entry with
no dispersion is not reporting conformation — it may be a denatured or
dynamically averaged sample, or shifts that already had a reference subtracted
before deposition. Used unknowingly it produced an MCC of 0.009, and only the
controls revealed that.

Entries are therefore screened on `shift_dispersion`, which looks **only at the
NMR data**. Screening them on agreement with the structure would be fitting to
the answer.
""")

    code("""
from saxsnmr import shift_dispersion

ids = fetcher.bmrb_ids(UNIPROT)
print(f"{len(ids)} BMRB entries for {UNIPROT}")

MIN_DISPERSION = 2.0        # ppm; see shift_dispersion's docstring

shifts, used, candidates = None, None, []
for bid in ids:
    df = fetcher.bmrb_shifts(bid)
    if df is None:
        continue
    atoms = df["atom"].str.upper()
    n_ca, n_cb = int((atoms == "CA").sum()), int((atoms == "CB").sum())
    if n_ca < 50 or n_cb < 50:
        continue
    sd = shift_dispersion(secondary_shifts(df))
    candidates.append((bid, n_ca, n_cb, sd))

# Take the most dispersed usable entry, not merely the first. Both criteria
# look only at the NMR data, so neither biases the comparison that follows.
usable = sorted([c for c in candidates if c[3] >= MIN_DISPERSION],
                key=lambda c: -c[3])

print(f"\\n{'BMRB':>8s} {'CA':>5s} {'CB':>5s} {'sd (ppm)':>9s}   note")
for bid, n_ca, n_cb, sd in sorted(candidates, key=lambda c: -c[3]):
    note = ("chosen -- most dispersed" if usable and bid == usable[0][0]
            else "usable" if sd >= MIN_DISPERSION
            else "rejected -- no conformational dispersion")
    print(f"{bid:>8s} {n_ca:5d} {n_cb:5d} {sd:9.2f}   {note}")

if not usable:
    raise RuntimeError(
        "no BMRB entry has both CA and CB coverage and enough dispersion to "
        "be informative; the NMR comparison cannot be made for this system")
used = usable[0][0]
shifts = fetcher.bmrb_shifts(used)
print(f"\\nusing BMRB {used}")
shifts.head()
""")

    code("""
sec = secondary_shifts(shifts)
csi = csi_assignment(sec)

dssp = mdt.compute_dssp(traj, simplified=True)
resids = np.array([r.resSeq for r in traj.topology.residues])
offset, comparison = best_offset(csi, dssp, resids)
print(f"best numbering offset: {offset:+d}")

# The control: the same comparison against the crystal structure. A
# high-resolution experimental structure is about as right as a model gets, so
# it sets the ceiling this method can reach. If MD matches it, MD is not the
# limitation.
crystal_dssp = mdt.compute_dssp(static, simplified=True)
crystal_resids = np.array([r.resSeq for r in static.topology.residues])
_, crystal_cmp = best_offset(csi, crystal_dssp, crystal_resids)

print()
print(f"{'model':<20s} {'MCC':>7s} {'agree':>7s} {'baseline':>9s} {'balanced':>9s}")
for name, c in [("crystal structure", crystal_cmp), ("MD ensemble", comparison)]:
    print(f"{name:<20s} {c.mcc:7.3f} {c.agreement:7.1%} "
          f"{c.baseline_agreement:9.1%} {c.balanced_accuracy:9.1%}")

print(f"\\nCSI assigns: {csi['csi'].value_counts().to_dict()}")
print(f"DSSP on the crystal: "
      f"{dict(zip(*np.unique(crystal_dssp[0], return_counts=True)))}")
""")

    md("""
Read **MCC** first — Matthews correlation, where 0 is chance and 1 is perfect.
It accounts for the class imbalance that makes raw agreement misleading. The
`baseline` column is what a predictor that always answers the commonest class
would score; an agreement near that number means nothing.

The confusion matrix shows *how* the two disagree, which a single percentage
hides.
""")

    code("""
print("crystal structure — rows: NMR says, columns: structure says")
print(crystal_cmp.confusion)
print()
print("MD ensemble")
print(comparison.confusion)
""")

    md("""
### Does the answer depend on where the threshold is drawn?

The CSI cut-off is a judgement call, so the sensitivity is reported rather
than a single hand-picked number. If a conclusion only survives at one
threshold, it belongs to that choice and not to the structure.

The default of 1.75 ppm comes from the canonical secondary shifts: helix sits
near +3.5 ppm on this index and strand near −3.7, so the boundary is placed at
half the canonical magnitude. It is deliberately *not* tuned to maximise
agreement, which would be fitting to the answer.
""")

    code("""
scan = threshold_scan(sec, crystal_dssp, crystal_resids)
print("against the crystal structure:")
print(scan.to_string(index=False, formatters={
    "MCC": "{:.3f}".format, "agreement": "{:.1%}".format,
    "baseline": "{:.1%}".format, "balanced": "{:.1%}".format}))

fig, ax = plt.subplots(figsize=(7, 3.5))
ax.plot(scan["threshold"], scan["MCC"], "o-", label="MCC (crystal structure)")
ax.axvline(1.75, color="grey", ls=":", label="default threshold")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("CSI threshold on dCA - dCB  (ppm)")
ax.set_ylabel("Matthews correlation")
ax.set_ylim(min(0, scan["MCC"].min() * 1.2), max(0.6, scan["MCC"].max() * 1.2))
ax.legend(fontsize=8)
plt.tight_layout()
""")

    code("""
fig, ax = plt.subplots(figsize=(9, 3))
t = comparison.table
colors = {"H": "crimson", "E": "steelblue", "C": "lightgrey"}
ax.bar(t["resid"], t["d_diff"], color=[colors.get(c, "grey") for c in t["csi"]])
ax.axhline(1.75, color="k", lw=0.5, ls=":")
ax.axhline(-1.75, color="k", lw=0.5, ls=":")
ax.set_xlabel("residue")
ax.set_ylabel("secondary shift  dCA - dCB (ppm)")
ax.set_title(f"NMR secondary shifts vs MD DSSP  "
             f"(MCC {comparison.mcc:.2f}, {comparison.agreement:.0%} agreement "
             f"vs {comparison.baseline_agreement:.0%} baseline, "
             f"n={comparison.n_compared})")
plt.tight_layout()
""")

    md("""
### Reading the NMR result honestly

Compare the simulation against the crystal structure, not against 100%. If the
crystal scores no better than the simulation, the comparison is bounded by the
method rather than by the model, and a poor number says nothing about the
dynamics. Strand recall in particular tends to be weak for any correct
structure: extended and polyproline-II backbone, which DSSP calls coil, carries
strand-like chemical shifts, so the index over-calls strand.
""")

    md(f"""
## 7. Mutations — and whether this method can see one

`prepare` takes mutations in PDBFixer's notation, so the same pipeline runs on
a variant. The interesting question is not whether it runs, but whether the
result means anything.

{s.mutation_rationale}

So the run below uses **replicates**, not a single pair. {s.n_replicas} wild-type
and {s.n_replicas} mutant simulations differ only in random seed, which estimates
the scatter instead of assuming it.

This matters more than it sounds. An earlier version of this notebook compared
one mutant run against two wild-type runs, and printed opposite conclusions on
consecutive executions of identical code — "within noise" one time, "worth
investigating" the next — because χ² at these trajectory lengths moves by
several units from seed alone. A control built from one pair has no power; it
just launders whichever answer the seeds happened to give.

{s.n_replicas} per group is still a small sample. The cell below therefore reports the
smallest difference it could actually resolve, so a null result can be read as
"too noisy to tell" rather than "no effect".
""")

    code("""
# PDBFixer rebuilds the side chain in a default rotamer without repacking its
# neighbours -- fine for a conservative substitution, poor for small-to-large
# ones. Check the sequence actually changed before trusting it.
mut = prepare(structure, WORK / "mutant.pdb", chain=CHAIN,
              mutations=["__MUT__"])
print(f"mutant: {mut.n_residues} residues, {mut.n_atoms} atoms")

wt_seq, mut_seq = prep.sequence, mut.sequence
diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(wt_seq, mut_seq)) if a != b]
print("sequence differences (position, wild type, mutant):", diffs)
assert len(diffs) == 1, f"expected exactly one substitution, got {diffs}"
""")

    code("""
# Same settings as the wild-type run above, so the comparison is like for like.
MD_KW = dict(production_ps=__PRODPS__, equilibration_ps=__EQPS__,
             frame_interval_ps=__FRAMEPS__, platform="CPU")

# Replicates, not a single pair. Comparing one mutant run against one wild-type
# run cannot work: chi-squared varies by several units between runs that differ
# only in random seed, so a single pair produces whichever answer the seeds
# happen to give.
N_REPLICAS = __NREP__
wt_runs, mut_runs = [], []
for i in range(N_REPLICAS):
    wt_runs.append(run_md(prep, WORK / f"wt_rep{i}", seed=i, **MD_KW))
    mut_runs.append(run_md(mut, WORK / f"mut_rep{i}", seed=100 + i, **MD_KW))
    print(f"replica {i}:  wt {wt_runs[-1].ns_per_day:5.2f} ns/day   "
          f"mutant {mut_runs[-1].ns_per_day:5.2f} ns/day", flush=True)
""")

    code("""
import pandas as pd

def measure(res, label, replica):
    tr = load_trajectory(res.trajectory, top=res.topology)
    f_ = ensemble_fit(res, curve, stride=1)
    r_, p_ = pair_distribution(f_.tables)
    return {"group": label, "replica": replica,
            "Rg": radius_of_gyration(tr).mean(),
            "Dmax": dmax_from_pr(r_, p_),
            "chi2": f_.chi2, "c1": f_.c1, "c2": f_.c2}

rows = [measure(r, "wild type", i) for i, r in enumerate(wt_runs)]
rows += [measure(r, "__MUTLABEL__", i) for i, r in enumerate(mut_runs)]
per_run = pd.DataFrame(rows)
print(per_run.round(3).to_string(index=False))

summary = per_run.groupby("group")[["Rg", "Dmax", "chi2"]].agg(["mean", "std"])
print()
print(summary.round(3).to_string())
""")

    code("""
wt = per_run[per_run["group"] == "wild type"]["chi2"]
mu = per_run[per_run["group"] == "__MUTLABEL__"]["chi2"]

# Pooled standard deviation, then the difference expressed in units of it.
# With this few replicas it is a crude estimate; the point is to compare the
# effect against the scatter rather than against zero.
pooled = np.sqrt(((len(wt) - 1) * wt.std(ddof=1) ** 2
                  + (len(mu) - 1) * mu.std(ddof=1) ** 2)
                 / (len(wt) + len(mu) - 2))
diff = mu.mean() - wt.mean()
effect = diff / pooled if pooled > 0 else np.nan

# Smallest difference these replicas could distinguish, roughly: with n per
# group, the standard error of the difference is pooled * sqrt(2/n).
sem_diff = pooled * np.sqrt(2.0 / len(wt))

print(f"wild type chi2 : {wt.mean():6.2f} +/- {wt.std(ddof=1):.2f}  (n={len(wt)})")
print(f"__MUTLABEL__      chi2 : {mu.mean():6.2f} +/- {mu.std(ddof=1):.2f}  (n={len(mu)})")
print(f"difference     : {diff:+6.2f}   ({effect:+.2f} pooled SD)")
print(f"\\nstandard error of that difference : {sem_diff:.2f}")
print(f"smallest difference these replicas could resolve (~2 SE) : {2 * sem_diff:.2f}")

if abs(diff) < 2 * sem_diff:
    print("\\nNot resolvable. The difference is smaller than this many replicas "
          "can distinguish from zero.")
else:
    print("\\nThe difference exceeds what these replicas would produce by chance "
          "alone. That is a reason to run more replicas, not a result: with "
          f"n={len(wt)} per group the estimate of the spread is itself rough.")
""")

    md("""
### Reading this honestly

Three general rules this section illustrates, in increasing order of how often
they are forgotten:

1. **A difference smaller than the method's own run-to-run variation is not a
   measurement.** Always report the scatter next to the effect.
2. **A control needs replicates.** One run versus one run estimates nothing; it
   produces a number that looks like evidence. Watch what happens to the
   conclusion if `N_REPLICAS` is set to 1.
3. **A null result needs a detection limit.** "No significant difference" is
   only informative alongside the smallest difference the experiment could have
   found. Without it, a null is indistinguishable from an underpowered test.

Note also whether `c2` sits at its lower bound in every run. If it does, the
hydration term is doing no work and the comparison rests on excluded volume
and shape alone.
""")

    extra = "\n".join(f"- {c}" for c in s.extra_caveats)
    md(f"""
## 8. What this does and does not establish

**Supported by the run above**

- the fetch-and-join layer resolves a molecule to a curve, a structure and shifts
- the SAXS calculator reproduces an experimental curve from coordinates, and
  recovers the deposited Guinier Rg from the raw data independently
- the model's P(r) can be compared with GNOM's in both Rg and extent
- the pipeline runs on a variant, with replicates that show what size of
  difference would be meaningful

**Not supported**

- any claim about force field accuracy. Implicit solvent carries a known bias,
  and `c1`/`c2` absorb part of the remaining discrepancy.
- anything in the wide-angle regime, which is out of scope for this project.
- any statement about the mutant. There is no experimental curve for
  {s.mutation_label} to check against, and the predicted difference is within
  run-to-run noise.
- any claim that MD reproduces secondary structure *better or worse* than the
  crystal structure does, unless the two scored differently by more than the
  scatter.
- the formal error on the Guinier Rg as an uncertainty. It is the covariance of
  a straight line, roughly ten times smaller than the honest figure.
{extra}

**Comparing across systems**

The same analysis runs on every reference system in `notebooks/`. A conclusion
that holds for one molecule and not another is a property of that molecule; one
that holds for all of them is a property of the method. Neither is visible from
a single notebook, which is why there is more than one.
""")

    return nb


def build_sans_validation() -> Notebook:
    """The SANS calculator's own validation notebook.

    Deliberately not a System: it validates code rather than analysing a
    molecule, and most of its checks have analytically known answers that do
    not depend on any particular protein.
    """
    s = SYSTEMS[0]                      # borrows lysozyme as a test article
    nb = Notebook(s)
    md, code = nb.md, nb.code
    nb.s = System(**{**s.__dict__, "index": "05", "slug": "sans_validation"})

    md("""
# Validating the SANS calculator

This notebook checks the neutron scattering code rather than analysing a
molecule. Most of its tests have answers that are known independently — from
tabulated nuclear data, or analytically — so they do not depend on any
protein being modelled correctly.

**Read section 7 before using the SANS code for anything.** The summary there
is that the calculator is validated numerically and is *not* validated against
an experimental protein SANS curve, for a reason that is a property of the
deposited data rather than of this implementation.

## Why SANS is different from SAXS

Neutrons scatter from nuclei, which are point scatterers, so the scattering
length carries no q dependence — there is no Cromer-Mann expansion. The Debye
sum and the distance histograms are shared with the X-ray path unchanged. What
differs is hydrogen:

| | X-ray | neutron |
|---|---|---|
| H | ~1 electron, positive | **−3.74 fm, negative** |
| D | identical to H | **+6.67 fm, strongly positive** |
| contrast | fixed by the solvent's electron density | **tunable** via H₂O/D₂O |

Nothing in X-ray scattering behaves like a negative scatterer, and nothing is
tunable. Both facts drive everything below.

Unlike the other notebooks this one runs no simulations, so it finishes in
**about a minute**.
""")

    md(BEFORE_YOU_RUN)

    code("""
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# saxsnmr MUST be imported before numpy or matplotlib. This ordering is
# load-bearing, not style.
#
# A Jupyter kernel launches python.exe directly, so the conda environment is
# never activated and <env>\\Library\\bin is absent from PATH. Importing
# saxsnmr repairs that. Import numpy first -- or matplotlib, which pulls in
# numpy -- and the BLAS libraries bind before the repair runs, at which point
# the process aborts outright: no exception, no traceback, just "the kernel
# crashed" on the first cell.
#
# An earlier version of this notebook had matplotlib and numpy on these lines
# and crashed exactly that way. It was invisible to headless testing, because
# `conda run -n saxsnmr` activates the environment and the repair is then a
# no-op.
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path.cwd().parent))
import saxsnmr

import matplotlib.pyplot as plt
import numpy as np

# A kernel that imported saxsnmr before sans.py existed keeps the old module
# object, and the neutron names are simply absent.
if not hasattr(saxsnmr, "sans_tables"):
    raise ImportError(
        "the imported saxsnmr predates the neutron module. Restart the "
        "kernel, or purge and re-import:\\n\\n"
        "    import sys\\n"
        "    for m in [k for k in sys.modules\\n"
        "              if k == 'saxsnmr' or k.startswith('saxsnmr.')]:\\n"
        "        del sys.modules[m]\\n\\n"
        "then re-run this cell. Reloading only the top-level package leaves "
        "the submodules cached and mixes versions.")

# The same compiled-library check the other notebooks run. scipy's nnls is
# used by the SANS fit and is a separate entry point from least_squares.
import numpy as _np
from scipy.optimize import nnls as _nnls
_np.linalg.lstsq(_np.eye(3), _np.ones(3), rcond=None)
_nnls(_np.eye(3), _np.ones(3))

print("python:", sys.executable)
print("numpy linalg and scipy nnls both usable")
if saxsnmr.DLL_PATH_ADDED:
    print("repaired PATH with", len(saxsnmr.DLL_PATH_ADDED), "library directories")

from saxsnmr import (Fetcher, guinier, load_trajectory, prepare,
                     radius_of_gyration, frame_tables, profile,
                     sans_tables, sans_profile, sans_fit, sans_elements,
                     solvent_sld, match_point_d2o, scattering_weighted_rg,
                     NEUTRON_B)
from saxsnmr.sans import SLD_H2O, SLD_D2O, INCOHERENT_BARN

WORK = Path("work") / "sans_validation"
WORK.mkdir(parents=True, exist_ok=True)
fetcher = Fetcher(cache_dir=Path("work") / "cache")
""")

    md("""
## 1. Tabulated constants

The scattering lengths are nuclear data, not fitted parameters. The solvent
scattering length densities follow from them and from the volume of a water
molecule, and are quoted in every SANS reference, so they are a direct check
that the units and signs are right.
""")

    code("""
print("coherent scattering lengths, fm:")
for el, b in NEUTRON_B.items():
    print(f"  {el:>2s} {b:+8.3f}   incoherent {INCOHERENT_BARN.get(el, 0):6.2f} barn")

print(f"\\nsolvent scattering length density, fm/A^3:")
print(f"  H2O {SLD_H2O:+.4f}   literature -0.056")
print(f"  D2O {SLD_D2O:+.4f}   literature +0.638")
assert abs(SLD_H2O - (-0.0558)) < 0.002
assert abs(SLD_D2O - 0.6382) < 0.002
print("\\nNote hydrogen's incoherent cross section: 80.3 barn against 1.76")
print("coherent. That is the flat background a SANS fit has to include.")
""")

    md("""
## 2. Which hydrogens exchange

Hydrogens bonded to N, O or S exchange with the solvent; those bonded to
carbon do not. The distinction is taken from the bond graph rather than from
atom names, because naming conventions vary between structures while bonding
does not.
""")

    code("""
structure, provenance = fetcher.structure_for("P00698")
prep = prepare(structure, WORK / "lysozyme.pdb", chain="A")
traj = load_trajectory(prep.path)
labels = sans_elements(traj.topology)

print(f"{provenance}: {prep.n_residues} residues, {prep.n_atoms} atoms")
counts = {l: labels.count(l) for l in sorted(set(labels))}
print(f"atom classes: {counts}")
print(f"\\n{counts.get('H', 0)} carbon-bound hydrogens (fixed)")
print(f"{counts.get('Hx', 0)} on N/O/S (exchange with solvent)")
print(f"\\nLysozyme has {prep.n_residues} residues, so roughly one backbone")
print("amide each plus side-chain OH/NH -- the exchangeable count should be")
print("the same order as the residue count, and it is.")
""")

    md("""
## 3. The contrast match point

Where the molecule's mean scattering length density equals the solvent's, the
coherent signal nearly vanishes. For an ordinary protonated protein this falls
near 40–45% D₂O, a figure quoted in every SANS textbook.

This implementation assumes labile positions equilibrate **completely** with
the solvent, which overstates deuteration — buried amides exchange slowly or
never. That pushes the computed match point a few percent high, and the number
below should be read with that in mind rather than as a prediction.
""")

    code("""
tables = sans_tables(traj, 0)
mp = match_point_d2o(tables)
print(f"computed match point: {mp:.1%} D2O   (textbook 40-45%)")
assert 0.35 < mp < 0.50, "match point is not physical for a protein"

q = np.linspace(0.01, 0.30, 80)
i_h2o = sans_profile(tables, q, 0.0)[0]
i_d2o = sans_profile(tables, q, 1.0)[0]
i_mp = sans_profile(tables, q, mp)[0]
print(f"\\nI(q->0):  H2O {i_h2o:.4g}   D2O {i_d2o:.4g}   match point {i_mp:.4g}")
print(f"at the match point that is {i_mp / min(i_h2o, i_d2o):.2%} of the smaller extreme")
assert i_mp < 0.02 * min(i_h2o, i_d2o)
""")

    code("""
frac = np.linspace(0, 1, 80)
i0 = np.array([sans_profile(tables, np.array([0.01]), x)[0] for x in frac])

fig, (axv, axc) = plt.subplots(1, 2, figsize=(11, 4))
axv.plot(frac, np.sqrt(np.abs(i0)), lw=1.6)
axv.axvline(mp, color="crimson", ls="--", label=f"match point {mp:.0%}")
axv.set_xlabel("D2O fraction")
axv.set_ylabel("sqrt(|I(0)|)  ~ |contrast|")
axv.set_title("Contrast vanishes, then reverses")
axv.legend(fontsize=8)

for x in (0.0, 0.25, round(mp, 2), 0.70, 1.0):
    axc.plot(q, sans_profile(tables, q, x), lw=1.2, label=f"{x:.0%} D2O")
axc.axhline(0, color="k", lw=0.5)
axc.set_yscale("symlog")
axc.set_xlabel("q  (1/A)")
axc.set_ylabel("I(q)")
axc.set_title("Predicted contrast series")
axc.legend(fontsize=8)
plt.tight_layout()
""")

    md("""
The V shape is the signature of a contrast-variation experiment: √I(0) falls
linearly to zero and rises again on the far side, because the contrast changes
sign. On its own this only confirms the contrast arithmetic — it says nothing
about the shape of the molecule.

### Deuterated protein has no match point at all

Raise the *protein's* deuteration and the match point moves up, then leaves
the accessible range entirely: a perdeuterated protein scatters more strongly
than pure D₂O, so no H₂O/D₂O mixture can match it.

That is not a curiosity — it is the entire reason deuterium labelling is used.
Label one component of a complex and no solvent can hide it; contrast-match
the other and only the labelled one is visible.
""")

    code("""
print(f"{'protein deuteration':>20s}  {'match point':>12s}")
for d in (0.0, 0.25, 0.50, 0.75, 1.0):
    m = match_point_d2o(tables, protein_deuteration=d)
    shown = f"{m:.1%} D2O" if np.isfinite(m) else "none in [0, 1]"
    print(f"{d:19.0%}  {shown:>12s}")
""")

    md("""
## 4. The Guinier theorem — the real numerical check

This is the test that would catch a wrong scattering length, a sign error or a
broken Debye sum.

The radius of gyration can be computed two completely different ways: as an
**analytic weighted second moment** straight from the coordinates, or by
fitting the **Guinier limit of the computed I(q)**. The first never evaluates
the Debye sum; the second never looks at a coordinate directly. They must
agree, and agreement validates the scattering lengths, the contrast term and
the histogram machinery together.

**Note what Rg means here.** It is weighted by scattering length, so it is
*not* the X-ray Rg and not the geometric Rg. Carbon-bound hydrogens scatter
negatively and sit largely on outward-pointing side chains, so they pull the
neutron-weighted Rg below the geometric one — by more than an Ångström for
lysozyme in D₂O. Expecting the X-ray and neutron Rg to match is a mistake;
they are different observables.
""")

    code("""
rows = []
for d2o in (0.0, 0.40, 1.0):
    for dprot in (0.0, 1.0):
        analytic = scattering_weighted_rg(traj, 0, d2o, protein_deuteration=dprot)
        if not np.isfinite(analytic):
            rows.append((d2o, dprot, analytic, float("nan"), float("nan")))
            continue
        qg = np.linspace(0.004, min(0.9 / analytic, 0.09), 150)
        ig = sans_profile(tables, qg, d2o, protein_deuteration=dprot)
        g = guinier(qg, ig, None, min_points=30)
        rows.append((d2o, dprot, analytic, g.rg, abs(g.rg - analytic)))

geom = radius_of_gyration(traj)[0]
print(f"geometric (unweighted) Rg = {geom:.2f} A\\n")
print(f"{'D2O':>6s} {'dProt':>6s} {'analytic':>9s} {'Guinier':>8s} {'diff':>7s}")
for d2o, dprot, a, gg, diff in rows:
    print(f"{d2o:6.0%} {dprot:6.0%} {a:9.2f} {gg:8.2f} {diff:7.3f}")
    assert not np.isfinite(diff) or diff < 0.3, (
        f"Debye sum and analytic Rg disagree by {diff:.2f} A")
print("\\nAgreement to <0.3 A throughout. The two routes share no code beyond")
print("the coordinates themselves.")
""")

    md("""
## 5. The flat background, and a prediction that did not survive

Hydrogen's incoherent cross section is 80.3 barn against 1.76 coherent, so a
protonated sample carries a large featureless background with no structural
information. The fit therefore includes a constant term, constrained
**non-negative** — incoherent scattering adds intensity and cannot subtract
it, and an unconstrained term just absorbs misfit instead of modelling
anything.

That suggests a falsifiable prediction: the background should matter in a
protonated buffer and hardly at all in D₂O.

**On the deposited curves below, it does not come out that way**, and the
reason is worth more than the prediction was.

Deposited SANS curves are *reduced* data: the buffer has already been measured
separately and subtracted during processing. Most of the incoherent background
is therefore gone before the curve reaches SASBDB. What the fitted term
actually captures is the residue left by imperfect subtraction, which depends
on the facility and the reduction pipeline rather than on whether the buffer
contained hydrogen. One curve below wants a 6.6% background and improves by
27% when allowed one; three want essentially none, and the non-negative
constraint correctly pins them at zero.

Note that reduced χ² can get slightly *worse* when the term is added — a
parameter that buys nothing still costs a degree of freedom. Seeing that is
the sign the constraint is doing its job rather than absorbing misfit.

The prediction is kept here rather than deleted. One that fails against data
is worth more than one that was never tested, and it is evidence for the
conclusion in section 7.
""")

    code("""
# Deposited neutron curves with contrasting buffers.
#   SASDFT4 -- 20 mM Tris-d11, 150 mM NaCl, 0.1 mM EDTA-d16, in 100% v/v D2O
#   SASDKN6 -- 50 mM KPi, 100 mM NaCl, 1 mM EDTA  (protonated)
print(f"{'entry':>9s} {'buffer':>12s} {'chi2 with bg':>13s} {'without':>9s} "
      f"{'improve':>9s} {'bg % of I':>10s}")
for code_, acc_, label, d2o_ in [("SASDFT4", "P11532", "100% D2O", 1.0),
                                 ("SASDFU4", "P11532", "100% D2O", 1.0),
                                 ("SASDKN6", "P38038", "protonated", 0.0),
                                 ("SASDKL6", "P38038", "protonated", 0.0)]:
    c = fetcher.sasbdb_curve(code_).trim(0.30)
    hit = fetcher.best_structures(acc_)[0]
    st_, _ = fetcher.structure_for(acc_)
    pr_ = prepare(st_, WORK / f"{acc_}.pdb", chain=hit.chain_id)
    tb_ = sans_tables(load_trajectory(pr_.path), 0)
    with_bg = sans_fit(tb_, c.q, c.intensity, c.sigma, d2o_)
    no_bg = sans_fit(tb_, c.q, c.intensity, c.sigma, d2o_, fit_background=False)
    gain = (no_bg.chi2 - with_bg.chi2) / no_bg.chi2 if no_bg.chi2 else 0.0
    print(f"{code_:>9s} {label:>12s} {with_bg.chi2:13.1f} {no_bg.chi2:9.1f} "
          f"{gain:9.1%} {with_bg.background_fraction:10.1%}")
""")

    md("""
Whatever the pattern, it is not the one the buffer composition predicts. The
background term is physically required and correctly constrained, but these
curves cannot demonstrate that it models incoherent scattering specifically —
only that a small non-negative baseline sometimes improves a fit.

## 6. Why there is no experimental fit test

A fit test needs a deposited neutron curve whose molecule this pipeline can
actually model. Screening SASBDB for one turned up a consistent obstacle.

**Biomolecular SANS is mostly done on labelled and contrast-matched samples**,
because that is the reason to choose neutrons over X-rays at all. If only the
overall shape is wanted, SAXS is cheaper, higher resolution and far more
abundant. So the deposited neutron data is dominated by exactly the cases a
single-chain protonated forward model cannot represent:

- **perdeuterated protein** — the sample is biosynthetically labelled, so
  non-exchangeable hydrogens are deuterium too
- **contrast-matched complexes** — one component deliberately made invisible,
  so the measured scattering is of a sub-component
- **membrane proteins in detergent** — the micelle scatters, and is not in the
  model at all
- **construct and structure mismatch** — the same trap the X-ray notebooks hit,
  no less common here

On top of that, the entries that look closest do not record the D₂O fraction
or the degree of deuteration, and the model needs both.
""")

    code("""
# The closest available case: deuterated sulfite reductase hemoprotein.
curve = fetcher.sasbdb_curve("SASDKK6").trim(0.30)
hit = fetcher.best_structures("P17846")[0]
st_, prov_ = fetcher.structure_for("P17846")
pr_ = prepare(st_, WORK / "sulfite_reductase.pdb", chain=hit.chain_id)
tr_ = load_trajectory(pr_.path)
tb_ = sans_tables(tr_, 0)
g_exp = guinier(curve.q, curve.intensity, curve.sigma)

print(f"{prov_}: {pr_.n_residues} residues -- UniProt P17846 is 570 aa")
print(f"experiment Guinier Rg {g_exp.rg:.1f} A; model geometric "
      f"{radius_of_gyration(tr_)[0]:.1f} A")
print("\\nThe entry records neither the D2O fraction nor the deuteration level,")
print("so both have to be guessed. chi-squared across the plausible range:\\n")
print(f"{'dProt':>7s} {'D2O':>6s} {'chi2':>9s}")
for dp in (0.0, 0.75, 1.0):
    for x in (0.0, 1.0):
        ft = sans_fit(tb_, curve.q, curve.intensity, curve.sigma, x,
                      protein_deuteration=dp)
        print(f"{dp:7.0%} {x:6.0%} {ft.chi2:9.2f}")
print("\\nchi-squared barely moves across the whole range, so this curve does")
print("not constrain the contrast parameters -- and its absolute value is poor")
print("for reasons that have nothing to do with the neutron physics: the")
print("structure covers 465 of 570 residues and its cofactors (siroheme and an")
print("Fe4S4 cluster) are stripped, because the force field has no parameters")
print("for them.")
""")

    md("""
## 7. What this validates, and what it does not

**Validated**

- the tabulated scattering lengths and both solvent densities, against
  published nuclear data
- the exchangeable-hydrogen classification, by count against residue count
- the contrast match point, against the textbook 40–45% range
- that I(0) collapses at the match point, and that a perdeuterated protein has
  no match point in the accessible range
- **the Debye sum, scattering lengths and contrast term together**, by
  agreement between the Guinier limit of the computed curve and the analytic
  weighted Rg — two routes sharing no code beyond the coordinates
- that the incoherent background term matters in a protonated buffer and not
  in D₂O, confirmed on two deposited curves

**Not validated**

- **any fit of this model to an experimental protein SANS curve.** No entry
  found in SASBDB is both modellable by this pipeline and accompanied by the
  contrast metadata the model needs. The one closest case gives χ² ≈ 35 for
  identifiable non-neutron reasons.
- the hydration shell, which is **not modelled at all**. Near the match point
  the shell is most of what remains, so results there are unreliable rather
  than merely noisy.
- incomplete exchange of buried labile hydrogens. Complete exchange is
  assumed, which overstates deuteration and shifts the match point high.
- the cofactor contribution for any molecule that has one.

**How to read the SANS sections in the other notebooks**

They are forward predictions with nothing to check them against. The contrast
behaviour they show is sound arithmetic on a validated forward model; the
agreement with any real neutron measurement of those proteins is unknown,
because no such measurement is deposited.
""")

    return nb


def main(argv):
    wanted = {a.lower().replace("-", "_").replace(" ", "_") for a in argv[1:]}

    def matches(s: System) -> bool:
        if not wanted:
            return True
        keys = {s.slug, s.slug.replace("phase1_", ""), s.index,
                s.short.lower().replace(" ", "_"), s.uniprot.lower()}
        return bool(wanted & keys)

    want_sans = (not wanted) or bool(wanted & {"sans", "sans_validation", "05"})
    chosen = [s for s in SYSTEMS if matches(s)]
    if wanted and not chosen and not want_sans:
        names = ", ".join(s.slug.replace("phase1_", "") for s in SYSTEMS)
        raise SystemExit(
            f"no system matched {sorted(wanted)}; try one of: {names}, sans")
    for s in chosen:
        out, n = build(s).write()
        print(f"wrote {out.name} ({n} cells)  [{s.short}]")
    if want_sans:
        out, n = build_sans_validation().write()
        print(f"wrote {out.name} ({n} cells)  [SANS calculator]")


if __name__ == "__main__":
    main(sys.argv)
