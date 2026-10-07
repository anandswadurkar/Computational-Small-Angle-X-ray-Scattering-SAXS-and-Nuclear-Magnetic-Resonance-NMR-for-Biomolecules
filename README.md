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
saxsnmr/md.py      PDB sanitising, structure prep, implicit-solvent MD
saxsnmr/nmr.py     secondary chemical shifts and comparison to DSSP
notebooks/02_phase1_lysozyme.ipynb   the pipeline end to end
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

Working: the fetch-and-join layer, structure sourcing and preparation, the SAXS calculator,
and the NMR comparison. Verified against lysozyme — the 0.65 A structure `2vb1` fitted to
SASBDB entry `SASDMJ2` gives reduced chi-squared **1.67** (c1 = 1.023, c2 = 0.007), with a
model Rg of 14.00 A against an experimental Guinier Rg of 13.98 A.

Not yet characterised: MD throughput on laptop-class hardware. Run the benchmark cell in
the notebook before planning trajectory lengths.

One caveat on the fitted parameters: `c1` and `c2` are defined by this implementation's own
parameterisation of excluded volume and hydration, so their values are not directly
comparable with the c1/c2 reported by CRYSOL or FoXS. Compare chi-squared between models
computed the same way, not against literature values from another program.
