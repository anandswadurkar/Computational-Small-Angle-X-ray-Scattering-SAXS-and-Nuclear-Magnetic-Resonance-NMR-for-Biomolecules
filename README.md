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
- NMR chemical shifts
- Short implicit-solvent MD, sized to run on a laptop
- Optional point mutations as exploratory input to the same pipeline, with a
  seed-to-seed control so a difference can be told from noise

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
  → shift predictor            → predicted vs assigned shifts
```

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
| implicit-solvent MD ensemble | **18.1 – 22.0** | 14.35 – 14.40 / 14.36 | 43.9 Å |

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
free to absorb anything: the MD ensemble reaches **43.9 Å** against an experimental Dmax
of 41.9 Å, while the crystal structure sits at 42.1 Å. The disagreement is not an artifact
of reciprocal-space fitting.

That is a real result, not a bug: the harness was built to detect this kind of
disagreement, and reporting it is the point. It does mean implicit solvent cannot support
a claim about force-field accuracy against SAXS.

### How large a χ² difference means anything

Running wild-type lysozyme twice, identical but for the random seed, gives χ² of 22.0 and
18.1 — a spread of **3.8** between two runs of the same molecule.

That is the noise floor for this configuration, and it is large. **Any χ² difference below
about 4 at 5 ps carries no information**, which is why the range above is quoted as a range
rather than a value. A 5 ps ensemble is not converged; it is long enough to show the
systematic expansion, which is several times this spread, and not much else.

This is also why the notebook's mutation section runs that control. E35Q — the catalytic
general acid, a charge change with almost no shape change — shifts χ² by 2.4, *less* than
the seed-to-seed spread. The correct conclusion is that the method cannot resolve it:

```
wild type (seed 0)   Rg 14.40 A   Dmax 43.9 A   chi2 22.0
wild type (seed 1)   Rg 14.35 A   Dmax 43.9 A   chi2 18.1
E35Q mutant          Rg 14.41 A   Dmax 44.1 A   chi2 19.6
```

There is also no experimental curve for E35Q in SASBDB, so there would be nothing to check
a prediction against even if one could be resolved. Simulating a variant is easy; knowing
whether the prediction is right is not.

NMR secondary-structure agreement is modest: 48.8% over 127 residues (helix 55.7%, strand
26.3%) against BMRB 4562. Strand agreement below chance partly reflects the simulation and
partly that the chemical shift index is weaker for strands than helices. These figures
move by a point or two between runs, for the same reason the χ² values do.

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
