# Computational Small Angle X-ray Scattering (SAXS) and Nuclear Magnetic Resonance (NMR) for Biomolecules

A locally deployable notebook pipeline that fetches experimental SAXS curves and NMR
chemical shifts for a protein of interest, runs short implicit-solvent molecular dynamics
from a matched structure, and reports how well the simulated observables reproduce the
measured ones.

Each notebook performs three tasks on one molecule:

1. **Short implicit-solvent MD**, producing an ensemble whose structural parameters can be
   compared against what a wet-lab instrument measures.
2. **Computational SAXS** — I(q) from coordinates, Guinier Rg, Kratky, and the distance
   distribution P(r) — checked against the deposited SASBDB curve and its GNOM transform.
3. **Computational NMR** — secondary chemical shifts converted to a per-residue secondary
   structure assignment, checked against BMRB, against the crystal structure, and against
   a no-information baseline.

This is a **validation harness**: it tests whether a given structure, force field and
workflow reproduce published experimental data. It is not a structure-determination or
structure-prediction tool. Its most useful output is a disagreement it can attribute.

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
- Small-angle **neutron** scattering, as a forward model: scattering lengths, H/D exchange,
  contrast variation and a non-negative incoherent background. Validated numerically but
  **not** against an experimental protein curve — see the SANS section below for why that
  is a property of the deposited data rather than of the code.

**Explicitly out of scope**

- Wide-angle scattering (WAXS). A uniform hydration shell is not valid near the water
  structure peak, and that regime needs explicit-solvent methods such as GROMACS-SWAXS
  or WAXSiS. Only ~3% of SASBDB entries reach q ≥ 1.0 Å⁻¹ anyway (see below).
- The **hydration shell in the neutron model**. The X-ray model fits a shell contrast;
  the neutron one does not model it at all, which matters most near the contrast match
  point where the shell is most of what remains.
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
  → neutron forward model      → contrast series and match point (prediction
                                 only; see the SANS section)
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

**The wanted chain is often not the first in the file.** `TER` and `SSBOND` records carry
their own chain identifiers, and filtering atoms without filtering those leaves a `TER`
for a discarded chain sitting before the first retained atom. OpenMM's PDB reader then
dies with `'NoneType' object has no attribute '_current_chain'`, which says nothing about
chains. This never fired on the reference structure, where the wanted chain happens to
come first — it only appeared when a second system needed chain `B`. A dangling `SSBOND`
to a removed chain is the same class of problem and is dropped too.

**GNOM writes P(r) in whatever length unit matches the input q.** A curve deposited in
1/nm has its r axis in nm, and the file carries no unit label. Assuming Ångström turns a
46 Å Dmax into 4.6 Å — small enough to look like a different molecule rather than an
obvious error. `Fetcher.pddf` converts using the same declared `angular_unit` the curve
uses, then **cross-checks the result against the entry's own reported Dmax and refuses
rather than guessing** if they disagree. One sampled entry fails that check, which is the
correct outcome: its deposited file and its metadata genuinely disagree.

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

## Reference systems

Three molecules run the identical analysis. A conclusion that holds for one and not the
others is a property of that molecule; one that holds for all three is a property of the
method. Neither is visible from a single notebook.

| | UniProt | SASBDB | structure | residues | crystal χ² | BMRB | D₂O match pt |
|---|---|---|---|---|---|---|---|
| **lysozyme** (hen egg white) | P00698 | SASDMJ2 | 2vb1, 0.65 Å | 129 | 1.67 | 39 | 46.8% |
| **RNase A** (bovine pancreatic) | P61823 | SASDPF7 | 6etk, 0.85 Å | 124 | 2.31 | 17 | 45.9% |
| **cytochrome c** (equine heart) | P00004 | SASDAB2 | 3jbt | 104 | 1.14 | 48 | 41.4% |

The computed neutron match points land at 41–47% D₂O across three unrelated proteins,
against a textbook 40–45%. Each is a separate calculation over a different composition, so
their clustering in the right band is a cross-system check the SANS model was not tuned to
pass. They are predictions: no deposited neutron curve exists for any of these three.

Lysozyme is the development and regression target: if the harness cannot reproduce
lysozyme, nothing else it reports is trustworthy. RNase A is the closest available
analogue — small, rigid, sub-Ångström structure — and exists to test whether the lysozyme
result is a property of lysozyme or of the method. Cytochrome c is a different fold class
and **carries a heme the force field cannot represent**, which makes it the most honest of
the three about what this pipeline cannot do.

### Finding these was harder than the 31% overlap suggests

Most obvious candidates fail, and they fail quietly:

- **α-synuclein** and **calmodulin**: dozens of SASBDB hits, but *zero* map to the
  accession. The join silently returns nothing.
- **Ubiquitin**: abundant data, but most entries are poly-ubiquitin chains — Rg 42–54 Å,
  Dmax 188–238 Å against a 12 Å monomer. The one monomeric entry found has no BMRB shifts.
- **SH3 domain**: the SIFTS "best structure" for the accession is a 318-residue complex,
  not the 62-residue domain. The crystal fit gives χ² ≈ 1000, which is the diagnostic
  working.
- **BPTI**: the deposited construct is 100 residues; the best structure resolves 41.
- Any **heme or metal protein**: `prepare` strips heterogens and amber14 has no parameters
  for them, so the simulation is silently of the apo chain.

Candidates were screened by running the whole chain — curve, Guinier, structure,
preparation, crystal fit — and keeping only what survived. A candidate that merely has
"SAXS data and NMR data" is not a candidate.

## Layout

```
saxsnmr/fetch.py   SASBDB, UniProt, PDBe/SIFTS, AlphaFold, BMRB, GNOM P(r)
saxsnmr/saxs.py    I(q) from coordinates, Guinier, P(r), chi-squared fitting
saxsnmr/sans.py    neutron scattering: scattering lengths, H/D exchange, contrast
saxsnmr/md.py      PDB sanitising, structure prep, implicit-solvent MD
saxsnmr/nmr.py     secondary chemical shifts and comparison to DSSP
tools/make_notebook.py                    generates the notebooks; edit this, not the .ipynb
notebooks/02_phase1_lysozyme.ipynb        reference system
notebooks/03_phase1_rnase_a.ipynb         same analysis, rigid analogue
notebooks/04_phase1_cytochrome_c.ipynb    same analysis, different fold + cofactor
notebooks/05_sans_validation.ipynb        validates the neutron calculator itself
```

The notebooks are **build output**. Editing cells in Jupyter works for a quick experiment,
but the generator overwrites them, so anything worth keeping belongs in the generator.
Hand-editing `.ipynb` JSON is error-prone and its diffs are unreadable, which is the whole
reason for generating it.

```
python tools/make_notebook.py             # all notebooks
python tools/make_notebook.py rnase_a     # just one system
python tools/make_notebook.py sans        # just the SANS validation
```

**The three molecule notebooks come from one `System` config each.** Only the molecule,
the accessions and the system-specific commentary differ — the analysis is literally the
same code. A harness that needed per-protein special-casing would not be a harness, so
where a system has its own trap (a cofactor, a chain that is not first in the file) it is
named in that system's config rather than worked around silently. The fourth notebook
validates the neutron calculator and is not a `System`, since it tests code rather than
analysing a molecule.

Adding another system means adding one `System` entry. Screen it end to end first; see
the reference-systems section for how many plausible candidates do not survive that.

### The generator refuses to write a broken notebook

Three checks run before any file is written. Each exists because the corresponding bug
shipped.

**Every code cell must compile.** Cell sources are built by string substitution into
templates, and a mangled escape produces a notebook that looks correct until someone runs
it — which happened, and put a `SyntaxError` in the first cell. Failing at generation
costs a second; failing in a notebook costs however long it took to reach that cell.

**`saxsnmr` must be imported before NumPy, SciPy, matplotlib, mdtraj or pandas.** The
generator names the offending line and refuses. See the Platform notes for why this
ordering is load-bearing; it crashed a kernel, and headless testing cannot detect it.

**The shared "Before you run this" note is emitted from one place** for all four
notebooks. The two builders drifting apart is exactly what produced the import-order bug:
the SANS notebook had been written with its own import cell that put NumPy first.

**The notebooks check that the imported `saxsnmr` is current.** A Jupyter kernel that
imported the package before it was last edited keeps the old module object, and new names
are simply absent — surfacing much later as an opaque `ImportError` in whichever cell
first needs one, after the simulations have already been paid for. The first cell now
verifies the API surface and, if it is stale, tells you to restart or to purge the cached
modules:

```python
import sys
for m in [k for k in sys.modules if k == "saxsnmr" or k.startswith("saxsnmr.")]:
    del sys.modules[m]
```

Purge the whole package, not just the top level: `importlib.reload(saxsnmr)` re-executes
`__init__.py` but leaves the submodules cached, so a new `sans.py` would call into a stale
`saxs.py` and fail more confusingly than before.

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

**Import `saxsnmr` before NumPy, SciPy, matplotlib, mdtraj or pandas.** This ordering is
load-bearing. A Jupyter kernel launches `python.exe` directly, so the conda environment is
never activated and `<env>\Library\bin` is absent from PATH; importing `saxsnmr` repairs
that. Import NumPy first -- or matplotlib, which pulls NumPy in -- and the BLAS libraries
bind before the repair runs, at which point the process aborts outright on the first cell:
no exception, no traceback.

This shipped in the SANS notebook and crashed a kernel. **Headless testing cannot catch
it**, because `conda run -n saxsnmr` activates the environment and makes the repair a
no-op; reproducing it requires invoking the environment's `python.exe` directly. The
generator now refuses to write a notebook whose first scientific import precedes
`saxsnmr`, and cell 1 of every notebook is tested under an unactivated interpreter.

**A crashing kernel is otherwise usually the machine, not the code.** When the SANS notebook killed
a kernel, measurement showed it peaks at 190 MB — the largest structure in it, 8,124 atoms,
costs 6 MB to prepare — while the machine had 1.5 GB free of 15.8 GB with the editor at
3.5 GB and a browser at 4.1 GB across 68 processes. Every operation suspected (`nnls` on
rank-deficient and all-zero inputs, `shrake_rupley` on the largest structure, a fit at the
match point where conditioning is worst) survived in an activated process. Measure the
notebook's footprint and the machine's free memory before suspecting the code; a long-lived
kernel still holding a previous notebook's state is the usual culprit, and restarting it
costs nothing.

**Every compiled entry point the notebooks use is checked in cell 1.** `scipy.optimize.nnls`
was missing from that check for a while: it is a separate compiled routine from
`least_squares`, used by the SANS fit, and anything that aborts there would otherwise do so
mid-notebook with no traceback.

**Almost every "broken library" symptom on Windows traces back to activation.** NumPy's
linear algebra, SciPy's optimisers and `mdtraj.load` all abort the process when the
environment is not activated, and all work normally when it is. Before suspecting a
package, confirm the environment is active.

## Status

All four notebooks run end to end, every code cell in one pass. The sections below are
long because this project records what it got wrong; the summary is here.

### At a glance

| | status |
|---|---|
| **fetch and join** | works — SASBDB, PDBe/SIFTS, RCSB, AlphaFold, BMRB, GNOM P(r) |
| **SAXS calculator** | validated: recovers deposited Guinier Rg to 0.01–0.09 Å on three systems |
| **implicit-solvent MD** | expands Rg by ~0.4 Å consistently; whether that helps **depends on the system** |
| **NMR comparison** | works, modest: MCC 0.39–0.44 against crystal structures, beating baseline |
| **mutation screening** | runs; **cannot resolve** any of the three mutations tested |
| **SANS calculator** | validated numerically; **no experimental protein fit** — see below |
| **tests, LICENSE, CITATION.cff** | **missing** — see Known gaps |

**The single most important result:** short implicit-solvent MD improves the SAXS fit for
RNase A (χ² 2.31 → 1.31), is mildly worse for cytochrome c, and much worse for lysozyme
(1.67 → 19.3). An earlier version of this README reported the lysozyme behaviour as a
property of the method. It is not.

**The most important caveat:** χ² is not comparable between these systems, because their
deposited uncertainties differ by 13–29×. Compare Rg in Ångström across systems; compare
χ² only within one.

### Known gaps

Measured against the ground rules in the parent `CLAUDE.md`:

- **No `tests/`.** The validation lives in notebooks and in throwaway scripts. The
  known-answer checks that exist (Guinier theorem, solvent SLDs, match point, lysozyme
  regression at χ² 1.662) should be `pytest` cases so a refactor cannot break them quietly.
- **No `LICENSE` or `CITATION.cff`.** Both are required by the ground rules; the licence
  choice is the author's.
- **No `pyproject.toml` or `pre-commit`.** `ruff` has never been run over this code.
- **No convergence study for the MD.** 5 ps is demonstrably not converged — the replicate
  scatter shows it — and no block averaging or autocorrelation-corrected error bars are
  computed.

### Lysozyme in detail

Lysozyme (`2vb1` against SASBDB `SASDMJ2`) is the development target:

| Model | reduced χ² | Rg (Guinier / P(r)) | Dmax from P(r) |
|---|---|---|---|
| experiment | — | 13.98 ± 0.26 deposited, **14.07 recomputed** / 14.06 | 41.9 Å |
| crystal structure `2vb1` | **1.67** | 14.00 / 13.91 | 42.1 Å |
| implicit-solvent MD ensemble | **19.3** | 14.36 / 14.31 | 43.6 Å |

### Across the three systems — the lysozyme result does not generalise

This is the finding that adding two more molecules produced, and it overturns what an
earlier version of this README claimed as a general conclusion.

| | experiment Rg | crystal Rg | MD Rg | crystal χ² | MD χ² |
|---|---|---|---|---|---|
| lysozyme | 14.06 Å | 13.91 | 14.31 | **1.67** | **19.3** |
| RNase A | 14.92 Å | 14.29 | **14.78** | **2.31** | **1.31** |
| cytochrome c | 13.00 Å | 12.85 | 13.24 | **1.14** | **1.44** |

(Rg from P(r), so no fitted parameter is involved.)

**For RNase A, implicit-solvent MD makes the fit better**, not worse — χ² 2.31 → 1.31. For
cytochrome c it is mildly worse. Only for lysozyme does it degrade badly. "Short
implicit-solvent MD makes SAXS agreement worse" was a statement about lysozyme that got
written up as a statement about the method.

**One mechanism explains all three.** GB implicit solvent expands the protein by about
0.4 Å in Rg — +0.40, +0.49, +0.39 for the three systems, which is strikingly consistent.
Whether that helps depends entirely on how far the crystal structure starts *below* the
solution value:

- RNase A's crystal is **0.63 Å too compact**, so the expansion moves it from 14.29 to
  14.78 against an experimental 14.92. The improvement is real but partly luck — the bias
  happened to point the right way and be roughly the right size.
- lysozyme's and cytochrome c's crystals are both only 0.15 Å too compact, so the same
  expansion overshoots and makes them worse.

The expansion is a consistent property of the solvent model. The *benefit* is a property
of the starting structure, which is why it cannot be generalised from one protein.

### χ² is not comparable between systems

| | points | median σ/I |
|---|---|---|
| lysozyme `SASDMJ2` | 147 | **0.0085** |
| RNase A `SASDPF7` | 1789 | 0.108 |
| cytochrome c `SASDAB2` | 1760 | 0.247 |

Lysozyme's deposited uncertainties are **13–29× tighter** than the others'. Since χ²
scales as 1/σ², the same relative model error costs lysozyme roughly 170× more χ² than
RNase A and 840× more than cytochrome c.

So lysozyme's χ² of 19.3 does **not** mean its model is thirteen times worse than
cytochrome c's at 1.44. It means lysozyme's measurement is far more precise, and the
harness is being held to a correspondingly higher standard. Rg in Ångström is the
quantity to compare across systems; χ² is only meaningful *within* one.

This also reframes the lysozyme result: a 0.25 Å overshoot in Rg, amplified by unusually
tight error bars, is what produced the dramatic-looking χ² of 19.3.

### How many replicates are enough? Not three, for lysozyme

The mutation test was run as three independent triplicates on lysozyme over the course of
this work. The verdicts were not stable:

| triplicate | wild-type χ² | difference | detection limit | printed verdict |
|---|---|---|---|---|
| 1 | 19.94 ± 1.54 | +0.41 | 3.57 | not resolvable |
| 2 | 18.43 ± 0.75 | +1.78 | 1.23 | **exceeds threshold** |
| 3 | 23.81 ± 5.63 | −3.73 | 7.73 | not resolvable |

The estimated scatter itself varies from 0.75 to 5.63 — a factor of seven. Two of three
triplicates say not resolvable and one says otherwise, which is what an underpowered test
looks like: the single positive verdict is consistent with chance, and should not be read
as an effect.

**n=3 is adequate for RNase A and cytochrome c and not for lysozyme**, because lysozyme's
tight error bars make its χ² both larger and more volatile. The notebook's wording handles
this correctly — when a difference exceeds the limit it says that is "a reason to run more
replicas, not a result" — but the honest summary is that this comparison needs more
replicates than the notebook runs by default.

Final run, all three systems:

| | wild type χ² | mutant χ² | difference | resolvable below | verdict |
|---|---|---|---|---|---|
| lysozyme | 23.81 ± 5.63 | 20.08 ± 3.63 (E35Q) | −3.73 | 7.73 | not resolvable |
| RNase A | 1.43 ± 0.14 | 1.41 ± 0.14 (H12A) | −0.02 | 0.22 | not resolvable |
| cytochrome c | 1.46 ± 0.16 | 1.35 ± 0.05 (M80A) | −0.11 | 0.20 | not resolvable |

All three mutations are chemically drastic and structurally conservative, and none is
resolvable. That is the expected and correct outcome.

### The data checks out; the simulation does not (lysozyme)

Recomputing the Guinier Rg from the raw curve gives 14.07 Å against the deposited
13.98 ± 0.26 Å — 0.4σ, over a self-selected window of 64 points spanning qRg 0.11–1.29.
The P(r) transform independently gives 14.06 Å. The experimental side of the comparison
is sound, so disagreement is the model's.

**For lysozyme specifically, short implicit-solvent MD makes agreement with SAXS
substantially worse.** The simulation expands the protein away from a crystal structure
that already matched experiment, and SAXS is most sensitive to exactly that. The expansion
appears within the first 2 ps and is unchanged across nonbonded cutoffs of 2.0, 1.5, 1.2
and 1.0 nm — and unchanged when the missing Debye screening is supplied (see below) — so it
is a property of the solvent model rather than a settings artifact. The fitted hydration
parameter `c2` falls to its lower bound, consistent with a model that is already too large.

See the cross-system section above before reading this as a general result: the same
expansion *improves* RNase A.

The P(r) comparison confirms this in real space, where no scale factor, `c1` or `c2` is
free to absorb anything: the MD ensemble reaches **43.6 Å** against an experimental Dmax
of 41.9 Å, while the crystal structure sits at 42.1 Å. The disagreement is not an artifact
of reciprocal-space fitting.

The expansion is also far more reproducible than χ² is. Across six independent runs Rg
stays within 14.34–14.43 Å (SD 0.03) while χ² ranges over 17.6–22.9. Rg is the better
quantity to reason about here; χ² inherits noise from every point in the curve.

That is a real result, not a bug: the harness was built to detect this kind of
disagreement, and reporting it is the point. It does mean implicit solvent cannot support
a claim about force-field accuracy against SAXS.

### Ionic strength: a real configuration error that is *not* the cause

`run_md` leaves Debye screening at OpenMM's default, which is **zero** — the simulation has
completely unscreened intramolecular electrostatics on a protein carrying net +8 e, while
`SASDMJ2` was measured in 150 mM NaCl at pH 3.8. That is a genuine mismatch with the
experiment and it was worth testing as an explanation.

It is not the explanation. Two seeds per condition:

| condition | Rg | net charge |
|---|---|---|
| experiment | **14.07 Å** | — |
| pH 7.0, 0 mM (the default) | 14.40 | +8 |
| pH 7.0, 150 mM screening | 14.34 | +8 |
| pH 3.8, 150 mM (matches the experiment) | **14.48** | +18 |

Supplying the screening moves Rg by 0.03 Å against a replicate scatter of ±0.025 — about
one standard deviation, so not a resolved effect. Matching the experimental pH makes
agreement *worse*, because lysozyme carries +18 e at pH 3.8. A model that fails harder as
the setup becomes more faithful is failing for physical reasons.

Two API notes, since this was harder to find than it should be: `implicitSolventSaltConc`
is **not** accepted on the `ForceField` path — it raises "specified but never used", which
reads like a typo rather than a wrong code path. `implicitSolventKappa` is the one that
works, and passing 0 reproduces the unset energy exactly, confirming the default really is
no screening.

### How the replicate control came to exist

Worth recording, because the failure was subtle and the fix was not obvious.

The first version estimated the scatter from **two** runs and reported a "noise floor" of
3.8. Re-running the identical code gave 1.75 for the same quantity — and the mutation
comparison built on top of it flipped from "within noise" to "worth investigating" between
consecutive executions. A control made of one run versus one run estimates nothing; it
converts whichever answer the seeds happened to give into something that reads like
evidence, in a section whose entire purpose was to prevent that.

Replicates and a stated detection limit replaced it. **"No significant difference" is
uninformative without the limit**, because a null result and an underpowered test look
identical; the limit turns a null into a bound on the effect rather than a claim about it.

The three-triplicate table above is what that machinery then revealed about its own
adequacy — which a single run would never have shown.

There is also no experimental curve for any of these mutants in SASBDB, so there would be
nothing to check a prediction against even if one could be resolved. Simulating a variant
is easy; knowing whether the prediction is right is not.

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

### A third defect: the BMRB entry has to be chosen, not taken

Running the analysis on two more proteins exposed one more problem. Selecting "the first
BMRB entry with more than 50 CA shifts" picked, for RNase A, an entry whose secondary
shifts have a standard deviation of **0.59 ppm** where five other entries for the same
protein give 3.6–5.0 ppm. An entry with no dispersion is not reporting conformation — it
may be a denatured or dynamically averaged sample, or shifts that already had a reference
subtracted before deposition.

The consequence was an MCC of **0.009**: no information at all. The controls caught it
(agreement 44% against a 98% majority baseline), which is exactly what they are for, but
the selection should not have offered it in the first place.

Entries are now screened on `shift_dispersion` and the **most dispersed** usable one is
taken. Both criteria look only at the NMR data; screening on agreement with the structure
would be fitting to the answer. The notebooks print the full candidate table with the
rejected entries visible.

The result is consistent across all three systems, which it was not before:

| | BMRB | MCC | agreement | baseline |
|---|---|---|---|---|
| lysozyme | 4562 | 0.440 | 63.8% | 47.2% |
| RNase A | 4031 | **0.401** (0.009 before) | 60.5% | 37.1% |
| cytochrome c | 17847 | **0.391** (0.242 before) | 63.5% | 46.2% |

An MCC near 0.4 is modest but real, and it is the ceiling this method reaches on a
correct structure. The MD ensembles score within ~0.02 of their crystal structures in
every case, so this comparison bounds the method and cannot separate the models.

### Neutron scattering: what the calculator is and is not

`saxsnmr/sans.py` implements SANS as a forward model. It is computationally *simpler* than
the X-ray case — nuclei are point scatterers, so the scattering length carries no q
dependence and there is no Cromer-Mann expansion — and it reuses the X-ray distance
histograms unchanged. What differs is hydrogen.

**Validated** (`notebooks/05_sans_validation.ipynb`):

| check | result |
|---|---|
| solvent SLDs vs published nuclear data | H₂O −0.0558, D₂O +0.6382 fm/Å³ — exact |
| exchangeable-H classification (lysozyme) | 696 carbon-bound, 263 on N/O/S |
| contrast match point | 46.8% D₂O (textbook 40–45%) |
| I(0) at the match point | 0.001% of either extreme |
| perdeuterated protein | **no match point exists** in [0, 1] |
| **Guinier limit vs analytic weighted Rg** | **agree to < 0.1 Å** across contrasts |

The last row is the real test. Rg can be computed either as an analytic weighted second
moment straight from coordinates, or by fitting the Guinier limit of the computed I(q).
The two share no code beyond the coordinates, so agreement validates the scattering
lengths, the contrast term and the Debye sum together.

**Rg is weighting-dependent, and this caught an error of mine.** I first asserted the
neutron and X-ray Guinier Rg should match and treated the 1.4 Å gap as a bug. It is not:
X-rays weight by electron count where hydrogen is ~1 and everything is positive, while
neutrons weight carbon-bound hydrogen at **−3.74 fm**. Those hydrogens sit largely on
outward-pointing side chains, so their negative contribution at large radius lowers the
neutron-weighted Rg. For lysozyme in D₂O it sits 1.25 Å below the geometric value. SANS
and SAXS Rg are different observables.

**Not validated: any fit to an experimental protein SANS curve.** This is the honest
limitation, and it is about the data.

Biomolecular SANS is mostly done on **labelled and contrast-matched** samples, because
that is the reason to choose neutrons at all — for overall shape, SAXS is cheaper, higher
resolution and far more abundant. So deposited neutron data is dominated by exactly the
cases a single-chain protonated forward model cannot represent. Screening SASBDB found:

- **zero** neutron entries for lysozyme or RNase A; the one "cytochrome c" hit is a
  three-component membrane complex in Fos-choline detergent
- perdeuterated protein (`SASDKK6`), contrast-matched heterodimers (`SASDKM6`, `SASDKN6`)
- membrane proteins in detergent micelles the model does not represent
- construct/structure mismatch, no less common than on the X-ray side

Perdeuteration support was added because of this — without it the model gets the sign of
the hydrogen contribution wrong for about a third of the atoms in a labelled sample. The
closest usable case (`SASDKK6`, deuterated sulfite reductase hemoprotein) still gives
χ² ≈ 35 for non-neutron reasons: the structure covers 465 of 570 residues and its siroheme
and Fe₄S₄ cluster are stripped. Its χ² also moves by under 3% across the whole plausible
range of D₂O fraction and deuteration, so it does not constrain the contrast parameters —
and the entry records neither.

**A prediction that failed.** I expected the incoherent background term to matter in a
protonated buffer and not in D₂O. It does not come out that way: of four curves, one
protonated case improves 27% and two are unaffected, while both D₂O cases improve
slightly. The explanation is that deposited curves are *reduced* data with the buffer
already subtracted, so what the term captures is residual imperfect subtraction — a
property of the reduction pipeline, not of the buffer. The failed prediction is kept in
the notebook; one that fails against data is worth more than one never tested.

The background is constrained **non-negative**, which is physics rather than convenience:
incoherent scattering adds intensity and cannot subtract it. An unconstrained term went
negative on three of five curves and was absorbing misfit rather than modelling anything.

**Two tools evaluated and not used**

- [NuMagSANS](https://github.com/AdamsMP92/NuMagSANS) computes SANS from **magnetic
  materials** — magnetization fields imported from MuMax3 or OOMMF, not atomic
  coordinates. Its own documentation states it does not handle biomolecular SANS. Right
  acronym, different physics.
- [Geant4](https://github.com/Geant4/geant4) simulates "the passage of particles through
  matter" for high-energy, nuclear, accelerator and medical physics. It is **not a
  molecular dynamics engine**: no force fields, no conformational dynamics, nothing that
  produces a protein ensemble. Geant4-DNA models radiation damage at nanoscale, which
  answers a different question from the one this project's MD asks.

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
