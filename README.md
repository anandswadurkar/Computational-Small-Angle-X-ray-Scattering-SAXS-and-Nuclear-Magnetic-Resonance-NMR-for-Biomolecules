# Computational Small Angle X-ray Scattering (SAXS) and Nuclear Magnetic Resonance (NMR) for Biomolecules

A locally deployable notebook pipeline that fetches experimental SAXS curves and NMR
chemical shifts for a protein of interest, runs short implicit-solvent molecular dynamics
from a matched structure, and reports how well the simulated observables reproduce the
measured ones.

This is a **validation harness**: it tests whether a given structure, force field and
workflow reproduce published experimental data. It is not a structure-determination or
structure-prediction tool.

## Scope

**In scope**

- Small-angle X-ray scattering, fitted over the q range where the data actually lives
- Standard SAXS diagnostics computed from the raw curve rather than taken on trust:
  Guinier with a self-consistent range, Kratky, and P(r) against the deposited GNOM
  transform
- NMR chemical shifts, scored against a majority-class baseline and the crystal
  structure so that an agreement figure can be interpreted
- Short implicit-solvent MD, sized to run on a laptop
- Optional point mutations as exploratory input to the same pipeline, run as replicates
  against a wild-type control so an effect can be told from seed noise, and reported with
  the smallest difference the replicates could actually resolve

**Explicitly out of scope**

- Wide-angle scattering (WAXS). A uniform hydration shell is not valid near the water
  structure peak, and that regime needs explicit-solvent methods such as GROMACS-SWAXS
  or WAXSiS. Only ~3% of SASBDB entries reach q ≥ 1.0 Å⁻¹ anyway (see below).
- Neutron scattering (SANS)
- Raw detector data, diffraction images, cryo-EM movies, or mass spectrometry

## Pipeline

```
molecule name
  → SASBDB search              → entry code
  → SASBDB entry summary       → UniProt + sequence + experimental I(q)
  → Guinier + Kratky           → Rg recomputed from the raw data; is it folded?
  → SASBDB GNOM file           → experimental P(r)
  → PDBe best_structures       → ranked experimental PDB structure
                                 (AlphaFold DB as fallback)
  → BMRB UniProt mapping       → assigned chemical shifts
  → PDBFixer                   → structure prep, optional mutation
  → OpenMM                     → short implicit-solvent MD
  → SAXS calculator            → per-frame I(q), ensemble average, χ² vs experiment
  → P(r) from coordinates      → real-space check, no free parameters
  → CSI vs DSSP                → secondary structure, against two controls:
                                 the crystal structure and a majority-class
                                 predictor
```

There is no shift predictor in that chain. SPARTA+, SHIFTX2 and UCBShift all lack a usable
Windows build, so the NMR comparison goes through secondary chemical shifts instead — a
weaker claim, and one that runs anywhere.

## Data sources

| Source | Supplies | Access |
|---|---|---|
| [SASBDB](https://www.sasbdb.org/) | experimental I(q), Rg, Dmax, sample metadata | REST, per-entry only |
| [PDBe / SIFTS](https://www.ebi.ac.uk/pdbe/) | UniProt → ranked experimental structures | REST |
| [RCSB PDB](https://www.rcsb.org/) | atomic coordinates | REST |
| [AlphaFold DB](https://alphafold.ebi.ac.uk/) | predicted structure fallback | REST, by UniProt |
| [BMRB](https://bmrb.io/) | assigned chemical shifts | REST, UniProt-mapped |
| [UniProt](https://www.uniprot.org/) | the join key between all of the above | REST |

## Measured characteristics of the sources

Sampled October 2026. Numbers are from random entry sampling, not full enumeration.

**SASBDB joinability** (n = 273 entries)

| Property | Fraction |
|---|---|
| has `uniprot_code` | 78% |
| has a `pdb_link` | 15% |
| has a sequence | 97% |

Because only 15% of entries link a structure, coordinates come from PDBe/SIFTS keyed on
UniProt rather than from SASBDB itself.

**SASBDB q-range** (n = 202 curves, units normalised via the `angular_unit` field)

| Statistic | q (Å⁻¹) |
|---|---|
| median q_max | 0.444 |
| 90th percentile q_max | 0.730 |
| fraction reaching q ≥ 1.0 | 3.0% |
| fraction reaching q ≥ 2.0 | 1.5% |

This is why the project is SAXS-only: the wide-angle regime has almost no deposited data,
and implicit-solvent hydration models are not reliable there in any case.

**SAXS and NMR overlap**: 31% of SASBDB UniProt accessions sampled also have BMRB entries.
The overlap is weighted toward well-studied model proteins, which suits a validation
harness.

## Pitfalls this code has to handle

**Angular units are not consistent.** About 60% of SASBDB entries declare `1/A` and 40%
declare `1/nm`. Always read the `angular_unit` field; never infer units from magnitude.
Guessing from magnitude silently inflates apparent q_max by 10× for a large minority of
entries.

**Lengths are in nanometres even when q is not.** SASBDB reports q in whatever
`angular_unit` declares, but always reports Rg and Dmax in nm. Mixing the two rescales
every length by ten — and 1.4 versus 14 is not obviously wrong at a glance. This package
converts both to Ångström on the way in and names the fields `guinier_rg_A`,
`pddf_rg_A`, `pddf_dmax_A` so that code written against the old names fails loudly
instead of being quietly wrong.

**A shared UniProt accession does not mean a shared molecule.** Entries for the same
accession differ in construct, tags, truncation, pH, temperature, oligomeric state and
ligands. Lysozyme (P00698) is a concrete example: deposited sequence lengths across its
SASBDB entries include 129, 131, 149 and 213 residues, and `SASDMJ4` reports
Rg ≈ 303 Å / Dmax ≈ 1200 Å where the monomeric entries report Rg ≈ 14 Å / Dmax ≈ 42 Å —
almost certainly a fibril sample. Construct metadata must be surfaced for the user to
check, not silently joined away.

**A Guinier fit must choose its own range.** The approximation holds only while qRg is
below roughly 1.3, which is circular because the valid range depends on the Rg being
measured. `guinier()` extends the window until the fitted Rg puts the last point past the
limit, and flags a curve where no window qualifies rather than returning a number. Run
against `SASDMJ4` it reports `[INVALID]` instead of a confident Rg of 303 Å.

**Caching text on Windows corrupted it.** Downloads already carry CRLF, and writing them
with `Path.write_text` translated the LF again, producing a blank line between every real
one when read back. Parsers that skip blanks survived; the GNOM P(r) parser, which treats
a blank line as the end of a block, silently returned a single data point. Cached text is
now written with `newline=""`. **Caches created before this fix should be deleted** — the
affected files are still on disk and still wrong.

**SAXS fitting has free parameters.** Excluded volume and hydration contrast are fitted by
most SAXS calculators. Report the fitted values alongside χ², or improvements will look
more meaningful than they are.

**Percentage agreement lies when classes are imbalanced.** Secondary structure is roughly
40% helix / 10% strand / 50% coil, so a model that always answers one class scores far
above zero. This project reported an NMR agreement figure for some time before a
majority-class baseline showed a no-information predictor beat it. Every classification
result here now carries its own baseline and a Matthews correlation; a percentage on its
own is not evidence.

**A control of one run against one run is not a control.** χ² at 5 ps moves by several
units from the random seed alone, so a single pair reports whichever answer the seeds
produced — this project published "within noise" and "worth investigating" from identical
code on consecutive executions before replicates were added. Any comparison between two
simulated conditions needs replicates on both sides and a stated detection limit; without
the limit, a null result cannot be distinguished from an underpowered test.

**Implicit solvent biases compaction.** GBSA-family models tend toward over-compact
structures, and Rg is precisely what SAXS is most sensitive to. Treat χ² from
implicit-solvent runs as a workflow check, not a force-field verdict.

## Reference system

Lysozyme (UniProt **P00698**) is the development and regression target:

- 26 SASBDB entries matching "lysozyme", with 14 confirmed mapped to P00698
- 39 BMRB entries carrying assigned chemical shifts
- PDB **2vb1** at 0.65 Å resolution

If the harness cannot reproduce lysozyme, nothing else it reports is trustworthy.

## Layout

```
saxsnmr/fetch.py   SASBDB, UniProt, PDBe/SIFTS, AlphaFold, BMRB, GNOM P(r)
saxsnmr/saxs.py    I(q) from coordinates, Guinier, P(r), chi-squared fitting
saxsnmr/md.py      PDB sanitising, structure prep, implicit-solvent MD
saxsnmr/nmr.py     secondary chemical shifts and comparison to DSSP
tools/make_notebook.py               generates the notebook; edit this, not the .ipynb
notebooks/02_phase1_lysozyme.ipynb   the pipeline end to end
```

The notebook is **build output**. Editing cells in Jupyter works for a quick experiment,
but `python tools/make_notebook.py` overwrites it, so anything worth keeping belongs in
the generator. Hand-editing `.ipynb` JSON is error-prone and its diffs are unreadable,
which is the whole reason for generating it.

## Install

```
conda env create -f environment.yml
conda activate saxsnmr
```

`environment.yml` pins `ipykernel` to 6.29 and `numpy` below 2.4, and takes everything
from conda-forge. Both pins are there for reasons documented in that file; relaxing them
reintroduces crashes that are hard to attribute. In particular, avoid `pip install`-ing a
compiled package into this environment on top of its conda build.

**On Windows, activation matters.** Conda keeps its compiled libraries in
`<env>\Library\bin`, which only joins PATH on activation. Invoking the environment's
`python.exe` by its full path instead makes NumPy's linear algebra abort the interpreter —
a hard process exit with no Python traceback, which is extremely misleading to diagnose.

Importing `saxsnmr` repairs this automatically, so scripts and notebooks work either way.
For anything that uses NumPy or OpenMM *before* importing `saxsnmr`, activate first with
`conda activate saxsnmr` or `conda run -n saxsnmr`.

Note that `os.add_dll_directory` is **not** a sufficient substitute: it does not cover the
transitive dependencies of the BLAS libraries. PATH is what works.

### Jupyter kernel

A Jupyter kernel launches `python.exe` directly, so it does not get the activated PATH.
That used to kill the kernel partway through the notebook with only "the kernel crashed"
and no traceback.

**Importing `saxsnmr` now repairs PATH itself**, so any kernel pointing at the
environment's interpreter works, including one made with
`python -m ipykernel install --user --name saxsnmr`. Import it before anything that
touches NumPy's linear algebra — the notebook's first cell does exactly that, and prints
how many directories it added.

`saxsnmr.DLL_PATH_ADDED` lists the repair for diagnostics; it is empty when the
environment was already activated.

No SAXS calculator binary is required. `saxs.py` implements the scattering model directly,
because neither FoXS nor Pepsi-SAXS offers a usable Windows build and CRYSOL needs a
licence that cannot be redistributed.

## Platform notes

**Use the CPU platform.** On the development machine (Intel Iris Xe integrated graphics),
`python -m openmm.testInstallation` reports OpenCL forces differing from CPU by 0.13, where
CPU and Reference agree to 6e-06. A platform computing wrong forces still produces a
trajectory that looks plausible. Verify this on your own hardware before using any platform
other than CPU.

**Almost every "broken library" symptom on Windows traces back to activation.** NumPy's
linear algebra, SciPy's optimisers and `mdtraj.load` all abort the process when the
environment is not activated, and all work normally when it is. Before suspecting a
package, confirm the environment is active.

## Status

The pipeline runs end to end — all 20 code cells of `02_phase1_lysozyme.ipynb` in one
pass. Measured on lysozyme (`2vb1` against SASBDB `SASDMJ2`):

| Model | reduced χ² | Rg (Guinier / P(r)) | Dmax from P(r) |
|---|---|---|---|
| experiment | — | 13.98 ± 0.26 deposited, **14.07 recomputed** / 14.06 | 41.9 Å |
| crystal structure `2vb1` | **1.66** | 14.00 / 13.90 | 42.1 Å |
| implicit-solvent MD ensemble | **19.9 ± 1.5** (n=3) | 14.34 – 14.39 / 14.36 | 43.4 Å |

### The data checks out; the simulation does not

Recomputing the Guinier Rg from the raw curve gives 14.07 Å against the deposited
13.98 ± 0.26 Å — 0.4σ, over a self-selected window of 64 points spanning qRg 0.11–1.29.
The P(r) transform independently gives 14.06 Å. The experimental side of the comparison
is sound, so disagreement is the model's.

**Short implicit-solvent MD makes agreement with SAXS substantially worse, not better.**
The simulation expands the protein away from a crystal structure that already matched
experiment, and SAXS is most sensitive to exactly that. The expansion appears within the
first 2 ps and is unchanged across nonbonded cutoffs of 2.0, 1.5, 1.2 and 1.0 nm, so it is
a systematic bias of the GB solvent model rather than a settings artifact. The fitted
hydration parameter `c2` falls to its lower bound, consistent with a model that is already
too large.

The P(r) comparison confirms this in real space, where no scale factor, `c1` or `c2` is
free to absorb anything: the MD ensemble reaches **43.4 Å** against an experimental Dmax
of 41.9 Å, while the crystal structure sits at 42.1 Å. The disagreement is not an artifact
of reciprocal-space fitting.

The expansion is also far more reproducible than χ² is. Across six independent runs Rg
stays within 14.34–14.43 Å (SD 0.03) while χ² ranges over 17.6–22.9. Rg is the better
quantity to reason about here; χ² inherits noise from every point in the curve.

That is a real result, not a bug: the harness was built to detect this kind of
disagreement, and reporting it is the point. It does mean implicit solvent cannot support
a claim about force-field accuracy against SAXS.

### How large a χ² difference means anything

Three wild-type replicates, identical but for the random seed, give χ² = **19.94 ± 1.54**.
A 5 ps ensemble is not converged; it is long enough to show the systematic expansion,
which is several times this scatter, and not much else.

**An earlier version of this section got this wrong, in a way worth recording.** It
estimated the scatter from *two* runs and reported a "noise floor" of 3.8. Re-running the
identical code gave 1.75 for the same quantity — and the mutation comparison built on top
of it flipped from "within noise" to "worth investigating" between consecutive executions.
A control made of one run versus one run estimates nothing; it converts whichever answer
the seeds happened to give into something that reads like evidence.

The notebook now runs three replicates per group and reports the detection limit:

```
wild type chi2 :  19.94 +/- 1.54  (n=3)
E35Q      chi2 :  20.35 +/- 2.68  (n=3)
difference     :  +0.41   (+0.19 pooled SD)
smallest difference these replicas could resolve (~2 SE) : 3.57
Not resolvable.
```

The difference is 0.41, against the 5.92 the two-run version reported for the same
comparison. E35Q is the catalytic general acid — a charge change with almost no shape
change — so "not resolvable" is the expected and correct answer.

Note the detection limit. **"No significant difference" is uninformative without it**,
because a null result and an underpowered test look identical. These replicates could not
have found anything smaller than ~3.6 in χ², so this result bounds the effect rather than
excluding it.

There is also no experimental curve for E35Q in SASBDB, so there would be nothing to check
a prediction against even if one could be resolved. Simulating a variant is easy; knowing
whether the prediction is right is not.

### NMR: what a control revealed

An earlier version of this README reported "48.8% agreement, strand 26.3%" and read it as
a modest result partly reflecting the simulation. Adding two controls showed that reading
was wrong, and that the metric was hiding the answer.

**Control 1 — the crystal structure.** Compare the same NMR data against DSSP on `2vb1`
instead of the trajectory. A 0.65 Å structure is as correct as a model gets, so it sets
the ceiling. It scored 49.6%, statistically indistinguishable from the simulation's 48.8%.
The simulation was never the limiting factor.

**Control 2 — a majority-class predictor.** Answering "helix" for every residue scored
**62.2%**, beating both. Any metric a no-information model wins is not reporting on
structure.

Two defects, found in that order:

1. **The CSI threshold was applied at the wrong scale.** `±0.7 ppm` is the textbook value
   for a *single* Δδ index; this code applies it to the difference `ΔδCA − ΔδCB`, which has
   roughly twice the dynamic range because helix raises CA *and* lowers CB. Against a
   measured standard deviation of 4.65 ppm, a ±0.7 band is a sliver: it assigned coil to 10
   residues out of 129 where DSSP finds 63. The threshold is now 1.75 ppm, derived from
   half the canonical helix/strand magnitude rather than tuned against a reference.

2. **Raw agreement was the wrong metric**, and it was the one doing the damage. Matthews
   correlation is 0.31 even at the old threshold — the index carried real signal the whole
   time, and percentage agreement under class imbalance was simply unable to show it.

After both fixes, against the crystal structure:

| | MCC | agreement | majority baseline |
|---|---|---|---|
| threshold 0.7 (old) | 0.317 | 49.6% | 62.2% ✗ |
| **threshold 1.75 (now)** | **0.440** | **63.8%** | 47.2% ✓ |

MCC stays between 0.32 and 0.46 across thresholds from 0.7 to 3.5, so the conclusion does
not depend on where the cut is drawn. The notebook reports that scan rather than a single
number.

Strand recall remains poor (23%) — but the crystal structure is barely better, so this is
a property of the index, not of the simulation. Extended and polyproline-II backbone that
DSSP calls coil carries strand-like shifts, so the index over-calls strand against *any*
correct structure. Without the crystal control that looked like a simulation defect.

`compare_to_dssp` now returns `mcc`, `baseline_agreement`, `balanced_accuracy` and a
confusion matrix, and `ShiftComparison.__repr__` flags a result that fails to beat its own
baseline.

### MD throughput

Roughly **6.5–8.5 ns/day** for a 1960-atom protein on a 4-core/8-thread laptop CPU, about
eleven times what the obvious settings give. The defaults in `run_md` encode this; see its
docstring for why each was chosen. The largest single factor is the nonbonded cutoff, which
at 2 nm enclosed every pair in the protein and so cost a great deal while achieving
nothing.

Run the notebook's benchmark cell on your own hardware before planning trajectory lengths.

One caveat on the fitted parameters: `c1` and `c2` are defined by this implementation's own
parameterisation of excluded volume and hydration, so their values are not directly
comparable with the c1/c2 reported by CRYSOL or FoXS. Compare chi-squared between models
computed the same way, not against literature values from another program.
