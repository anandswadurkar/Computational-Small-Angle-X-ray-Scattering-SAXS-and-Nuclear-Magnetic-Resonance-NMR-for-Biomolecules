"""Generate the Phase 1 notebook.

The notebook is build output; this file is the source. Editing cells in
Jupyter and saving works for a quick experiment, but the next run of this
script overwrites it, so anything worth keeping belongs here.

    python tools/make_notebook.py

Hand-editing .ipynb JSON is error-prone and its diffs are unreadable, which is
the whole reason for generating it.
"""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

cells = []
_counter = [0]


def _next_id():
    _counter[0] += 1
    return f"cell{_counter[0]:03d}"


def md(text):
    cells.append({"cell_type": "markdown", "id": _next_id(), "metadata": {},
                  "source": text.strip("\n").splitlines(keepends=True)})


def code(text):
    cells.append({"cell_type": "code", "id": _next_id(), "execution_count": None,
                  "metadata": {}, "outputs": [],
                  "source": text.strip("\n").splitlines(keepends=True)})


md("""
# Phase 1 — validating short implicit-solvent MD against SAXS and NMR

Reference system: **hen egg white lysozyme**, UniProt `P00698`.

This notebook runs the full loop end to end:

1. survey the candidate SASBDB entries and **check they describe the same molecule**
2. load one experimental scattering curve, normalised to 1/A
3. **check the curve against itself** — recompute Guinier Rg from the raw data,
   and look at the Kratky plot before assuming the particle is folded
4. source atomic coordinates through UniProt
5. prepare and simulate in implicit solvent
6. compute an ensemble SAXS profile and fit it
7. compare the distance distribution P(r) with the deposited GNOM transform
8. compare NMR chemical shifts to the simulated secondary structure
9. run a point mutant **with a seed-to-seed control**, to see whether the
   method can resolve the change at all

Lysozyme is the regression target because it is exhaustively characterised.
If the harness cannot reproduce lysozyme, nothing else it reports is
trustworthy.
""")

md("""
> **Use a kernel from the `saxsnmr` environment.**
>
> On Windows, conda keeps its compiled libraries in `<env>\Library\bin`, which only
> joins the search path when the environment is activated. A Jupyter kernel launches
> `python.exe` directly, so it normally misses them, and the first call into NumPy's
> linear algebra then kills the kernel with no exception and no traceback — reported
> only as *"the kernel crashed"*.
>
> Importing `saxsnmr` repairs PATH itself, so the notebook works whichever way the
> kernel was started. That is why the cell below imports it before anything else.
""")

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

print("python:", sys.executable)
if saxsnmr.DLL_PATH_ADDED:
    print("repaired PATH with", len(saxsnmr.DLL_PATH_ADDED), "library directories")
else:
    print("PATH already complete (environment was activated)")

checks = {
    "numpy linalg":   "import numpy as np; np.linalg.lstsq(np.eye(3), np.ones(3), rcond=None)",
    "scipy optimize": "from scipy.optimize import least_squares; least_squares(lambda p: p - 3.0, [0.0])",
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
    print("If the kernel dies instead of printing this, it was started without")
    print("the environment's library directory on PATH -- see the README.")
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
from saxsnmr.nmr import secondary_shifts, csi_assignment, best_offset

WORK = Path("work")
WORK.mkdir(exist_ok=True)
fetcher = Fetcher(cache_dir=WORK / "cache")
UNIPROT = "P00698"
""")

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
codes = fetcher.search_sasbdb("lysozyme")
print(f"{len(codes)} candidate entries")

survey = summarise_constructs(fetcher, codes[:12])
survey[["code", "uniprot", "length", "guinier_rg_A", "pddf_dmax_A",
        "angular_unit", "type_of_curve"]]
""")

md("""
### Reading the survey

Monomeric lysozyme is 129 residues with Rg ~14 A and Dmax ~42 A. Entries
reporting a much larger Rg or Dmax are aggregates or fibrils rather than the
folded monomer, and must not be averaged together with it.
""")

code("""
CODE = "SASDMJ2"   # monomeric lysozyme

curve = fetcher.sasbdb_curve(CODE)
print(curve)
print(f"SASBDB Guinier Rg = {curve.guinier_rg_A:.2f} +/- {curve.guinier_rg_err_A:.2f} A")
print(f"SASBDB P(r):  Rg = {curve.pddf_rg_A:.2f} A, Dmax = {curve.pddf_dmax_A:.1f} A")

# Restrict to where an implicit hydration shell is valid (see README).
curve = curve.trim(q_max=0.50)
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
print(f"  SASBDB reports {curve.guinier_rg_A:.2f} +/- {curve.guinier_rg_err_A:.2f} A")
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
Lysozyme should show that clearly. If a future target does not, the premise
that one short trajectory represents its solution ensemble is already in
doubt, and a poor chi-squared later would be a symptom rather than a cause.

## 2. Coordinates

SASBDB links a structure for only about 15% of entries, so coordinates come
from PDBe/SIFTS keyed on UniProt, falling back to AlphaFold.

Preparation strips pre-existing hydrogens, resolves alternate conformations and
selects a single chain. Lysozyme's best structure, `2vb1` at 0.65 A, needs all
three: it resolves its own hydrogens, carries altloc codes A/B/C, and labels two
different chains `A`.
""")

code("""
for hit in fetcher.best_structures(UNIPROT)[:5]:
    print(f"  {hit.pdb_id} chain {hit.chain_id}  {hit.resolution} A  {hit.method}")

structure, provenance = fetcher.structure_for(UNIPROT)
print("\\nusing:", provenance)

prep = prepare(structure, WORK / "prepared.pdb", chain="A")
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

Implicit solvent is what puts this on a laptop. It also biases towards
over-compact structures, and Rg is exactly what SAXS measures most sensitively —
so treat the fit below as a check that the workflow is sound, not as a verdict
on the force field.
""")

code("""
# run_md's tuned defaults give about 6.5-8.5 ns/day for this system on a
# laptop CPU -- roughly eleven times the obvious settings, mostly from a
# shorter nonbonded cutoff. 5 ps therefore takes about a minute; 50 ps is
# around ten. Raise production_ps once the benchmark above tells you your rate.
result = run_md(prep, WORK, production_ps=5.0, equilibration_ps=2.0,
                frame_interval_ps=0.5, platform="CPU")
print(result)
""")

code("""
import mdtraj as mdt

# load_trajectory pairs the DCD with its topology and converts units.
traj = load_trajectory(result.trajectory, top=result.topology)
rg = radius_of_gyration(traj)

fig, ax = plt.subplots(figsize=(6, 3))
ax.plot(np.arange(len(rg)) * 0.5, rg, lw=1)
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

## 6. NMR: chemical shifts against simulated secondary structure

Predicting raw shifts needs SPARTA+, SHIFTX2 or UCBShift, none of which builds
on Windows. The portable route is secondary chemical shifts: subtracting
random-coil values leaves a signal whose sign reports helix versus strand, which
can be compared with DSSP over the trajectory.

This checks that the ensemble has the right secondary structure where NMR says
it should. It is a weaker claim than reproducing every shift, and it is the
claim this notebook can actually support.
""")

code("""
ids = fetcher.bmrb_ids(UNIPROT)
print(f"{len(ids)} BMRB entries for {UNIPROT}: {ids[:10]}")

shifts = None
for bid in ids:
    df = fetcher.bmrb_shifts(bid)
    if df is not None and {"CA", "CB"} & set(df["atom"].str.upper()):
        n_ca = (df["atom"].str.upper() == "CA").sum()
        if n_ca > 50:
            shifts, used = df, bid
            break
print(f"using BMRB {used}: {len(shifts)} assigned shifts, {n_ca} CA")
shifts.head()
""")

code("""
sec = secondary_shifts(shifts)
csi = csi_assignment(sec)
print(csi["csi"].value_counts().to_dict())

dssp = mdt.compute_dssp(traj, simplified=True)
resids = np.array([r.resSeq for r in traj.topology.residues])

offset, comparison = best_offset(csi, dssp, resids)
print(f"\\nbest numbering offset: {offset:+d}")
print(comparison)
""")

code("""
fig, ax = plt.subplots(figsize=(9, 3))
t = comparison.table
colors = {"H": "crimson", "E": "steelblue", "C": "lightgrey"}
ax.bar(t["resid"], t["d_diff"], color=[colors.get(c, "grey") for c in t["csi"]])
ax.axhline(0.7, color="k", lw=0.5, ls=":")
ax.axhline(-0.7, color="k", lw=0.5, ls=":")
ax.set_xlabel("residue")
ax.set_ylabel("secondary shift  dCA - dCB (ppm)")
ax.set_title(f"NMR secondary shifts vs MD DSSP "
             f"({comparison.agreement:.0%} agreement over {comparison.n_compared} residues)")
plt.tight_layout()
""")

md("""
## 7. Mutations — and whether this method can see one

`prepare` takes mutations in PDBFixer's notation, so the same pipeline runs on
a variant. The interesting question is not whether it runs, but whether the
result means anything.

Glu35 is lysozyme's catalytic general acid. **E35Q** is a classic inactivating
substitution: it removes the charge while keeping almost exactly the same
shape. That makes it a good test of the method's limits, because a change that
drastic for the enzyme's *function* should be nearly invisible to SAXS, which
only sees the shape.

So the run below has a control. Wild type is simulated a second time with a
different random seed, giving a baseline for how much two runs of the *same*
molecule differ. A mutant-versus-wild-type difference only means something if
it is larger than that.

Without this control, any difference at all looks like a result.
""")

code("""
# PDBFixer rebuilds the side chain in a default rotamer without repacking its
# neighbours -- fine for a conservative substitution like Glu -> Gln, poor for
# small-to-large ones. Check the sequence actually changed before trusting it.
mut = prepare(structure, WORK / "mutant.pdb", chain="A",
              mutations=["GLU-35-GLN"])
print(f"mutant: {mut.n_residues} residues, {mut.n_atoms} atoms")

wt_seq, mut_seq = prep.sequence, mut.sequence
diffs = [(i + 1, a, b) for i, (a, b) in enumerate(zip(wt_seq, mut_seq)) if a != b]
print("sequence differences (position, wild type, mutant):", diffs)
assert diffs == [(35, "E", "Q")], "the mutation did not apply as intended"
""")

code("""
# Same settings as the wild-type run above, so the comparison is like for like.
MD_KW = dict(production_ps=5.0, equilibration_ps=2.0,
             frame_interval_ps=0.5, platform="CPU")

mut_result = run_md(mut, WORK / "mutant_md", **MD_KW)
print("mutant      ", mut_result)

# The control: wild type again, identical except for the random seed.
wt2_result = run_md(prep, WORK / "wt_repeat", seed=1, **MD_KW)
print("wt (seed 1) ", wt2_result)
""")

code("""
rows = []
for label, res in [("wild type (seed 0)", result),
                   ("wild type (seed 1)", wt2_result),
                   ("E35Q mutant", mut_result)]:
    tr = load_trajectory(res.trajectory, top=res.topology)
    f_ = ensemble_fit(res, curve, stride=1)
    r_, p_ = pair_distribution(f_.tables)
    rows.append({"model": label,
                 "Rg (A)": radius_of_gyration(tr).mean(),
                 "Dmax (A)": dmax_from_pr(r_, p_),
                 "chi2": f_.chi2, "c1": f_.c1, "c2": f_.c2})

import pandas as pd
table = pd.DataFrame(rows).set_index("model").round(3)
print(table.to_string())

seed_spread = abs(table.loc["wild type (seed 1)", "chi2"]
                  - table.loc["wild type (seed 0)", "chi2"])
mutant_shift = abs(table.loc["E35Q mutant", "chi2"]
                   - table.loc["wild type (seed 0)", "chi2"])
print(f"\\nchi2 difference between two wild-type runs : {seed_spread:.3f}")
print(f"chi2 difference, mutant vs wild type        : {mutant_shift:.3f}")
print("\\nVERDICT:", "mutant effect is within run-to-run noise -- not resolvable"
      if mutant_shift <= seed_spread else
      "mutant effect exceeds the seed-to-seed baseline -- worth investigating")
""")

md("""
### Reading this honestly

With 5 ps and two free parameters, the expected outcome is that the mutant is
**indistinguishable** from the wild type. That is the correct result, not a
failure: E35Q is a charge change, SAXS at this resolution measures shape, and
two runs of the same molecule already differ by a comparable amount.

The general rule this illustrates: **a difference smaller than the method's own
run-to-run variation is not a measurement.** Report the baseline whenever
reporting a comparison, and note that the baseline here comes from a single
pair of short runs, which is itself a weak estimate of the spread.

There is also no experimental curve for E35Q in SASBDB, so there is nothing to
validate the mutant prediction against. Simulating a variant is easy; knowing
whether the prediction is right is not.

## 8. What this does and does not establish

**Supported by the run above**

- the fetch-and-join layer resolves a molecule to a curve, a structure and shifts
- the SAXS calculator reproduces an experimental curve from coordinates, and
  recovers the deposited Guinier Rg from the raw data independently
- the model's P(r) agrees with GNOM's in both Rg and extent
- the pipeline runs on a variant, with a control that shows what size of
  difference would be meaningful

**Not supported**

- any claim about force field accuracy. Implicit solvent biases compaction, and
  `c1`/`c2` absorb part of the remaining discrepancy.
- anything in the wide-angle regime, which is out of scope for this project.
- any statement about a mutant. There is no experimental curve for E35Q to
  check against, and the predicted difference is within run-to-run noise.
- the formal error on the Guinier Rg as an uncertainty. It is the covariance of
  a straight line, roughly ten times smaller than the honest figure.

**The headline negative result**

Short implicit-solvent MD makes agreement with SAXS *worse*, not better: the
crystal structure alone fits better than the ensemble derived from it. The
harness was built to detect that kind of disagreement, so reporting it is the
point rather than an embarrassment. See the README's Status section.
""")

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "Python (saxsnmr)", "language": "python",
                       "name": "saxsnmr"},
        "language_info": {"name": "python", "version": "3.11"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

out = REPO / "notebooks" / "02_phase1_lysozyme.ipynb"
out.write_text(json.dumps(nb, indent=1), encoding="utf-8")
print("wrote", out, f"({len(cells)} cells)")
