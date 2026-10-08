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
- NMR chemical shifts
- Short implicit-solvent MD, sized to run on a laptop
- Explicit-solvent MD on a GPU, added to answer a question the implicit-solvent result
  raised rather than planned from the start (see Status)
- Optional point mutations as exploratory input to the same pipeline

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
  → PDBe best_structures       → ranked experimental PDB structure
                                 (AlphaFold DB as fallback)
  → BMRB UniProt mapping       → assigned chemical shifts
  → PDBFixer                   → structure prep, optional mutation
  → OpenMM                     → short implicit-solvent MD
  → SAXS calculator            → per-frame I(q), ensemble average, χ² vs experiment
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

**A shared UniProt accession does not mean a shared molecule.** Entries for the same
accession differ in construct, tags, truncation, pH, temperature, oligomeric state and
ligands. Lysozyme (P00698) is a concrete example: deposited sequence lengths across its
SASBDB entries include 129, 131, 149 and 213 residues, and `SASDMJ4` reports
Rg ≈ 30 / Dmax ≈ 120 where the monomeric entries report Rg ≈ 1.4 / Dmax ≈ 4.2 — almost
certainly a fibril sample. Construct metadata must be surfaced for the user to check, not
silently joined away.

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
saxsnmr/fetch.py   SASBDB, UniProt, PDBe/SIFTS, AlphaFold, BMRB
saxsnmr/saxs.py    SAXS profile from coordinates, and chi-squared fitting
saxsnmr/md.py      PDB sanitising, structure prep, implicit and explicit solvent MD
saxsnmr/nmr.py     secondary chemical shifts and comparison to DSSP
notebooks/02_phase1_lysozyme.ipynb         the pipeline end to end, on a laptop
notebooks/03_colab_explicit_solvent.ipynb  explicit solvent on a Colab GPU
```

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

The pipeline runs end to end. Measured on lysozyme (`2vb1` against SASBDB `SASDMJ2`):

| Model | reduced chi-squared | Rg |
|---|---|---|
| experiment | — | 13.98 A (Guinier) |
| crystal structure `2vb1` | **1.67** | 14.00 A |
| implicit-solvent MD ensemble | **14.5 – 20.2** | 14.2 – 14.4 A |

**Short implicit-solvent MD makes agreement with SAXS substantially worse, not better.**
The simulation expands the protein away from a crystal structure that already matched
experiment, and SAXS is most sensitive to exactly that. The expansion appears within the
first 2 ps and is unchanged across nonbonded cutoffs of 2.0, 1.5, 1.2 and 1.0 nm, so it is
a systematic bias of the GB solvent model rather than a settings artifact. The fitted
hydration parameter `c2` falls to its lower bound, consistent with a model that is already
too large.

That is a real result, not a bug: the harness was built to detect this kind of
disagreement, and reporting it is the point. It does mean implicit solvent cannot support
a claim about force-field accuracy against SAXS. Explicit solvent, on hardware that can
afford it, is the route to that.

### Explicit solvent

`solvate` and `run_md_explicit` provide the explicit-solvent path — PME, NPT with a Monte
Carlo barostat, and a protein-only trajectory so water never reaches disk. Solvating
lysozyme takes it from 1,960 to 22,508 atoms in a 6.2 nm box, an 11.5x growth, which is why
this needs a GPU; `notebooks/03_colab_explicit_solvent.ipynb` runs it on a free Colab T4.

**No explicit-solvent result is reported here yet.** The code path is validated on CPU —
PME setup, the barostat, the protein-only subset reporter, checkpoint write, resume, and
that the output feeds `ensemble_fit` — but the CUDA platform and the physics have not been
exercised, because this machine has no usable GPU. Treat the numbers in that notebook's
comparison table as Phase 1's, with the explicit-solvent row still empty.

Two details are worth knowing before running it. Solvation is slow and single-threaded (48
minutes on the development laptop), so `solvate` reuses an existing box by default rather
than rebuilding it after a reconnect. And `production_ns` is a cumulative target: re-running
the MD cell resumes and extends rather than restarting, because Colab sessions drop.

NMR secondary-structure agreement is modest: 49.6% over 127 residues (helix 58.2%, strand
26.3%) against BMRB 4562. Strand agreement below chance partly reflects the simulation and
partly that the chemical shift index is weaker for strands than helices.

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
