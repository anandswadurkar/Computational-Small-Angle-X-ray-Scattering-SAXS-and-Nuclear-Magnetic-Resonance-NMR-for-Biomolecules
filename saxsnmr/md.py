"""Structure preparation and molecular dynamics, implicit or explicit solvent.

Implicit solvent (``run_md``) is what makes this project run on a laptop:
dropping explicit water takes a small protein from tens of thousands of atoms
to a couple of thousand.

The cost is a measured bias, not a hypothetical one. For lysozyme the GB model
expands the protein to a radius of gyration of about 14.2 A against an
experimental Guinier value of 13.98 A, reaching it within the first two
picoseconds and holding there across every nonbonded cutoff tried. Since SAXS
is most sensitive to exactly that quantity, chi-squared rises from 1.67 for the
crystal structure to 14-20 for the simulated ensemble. Treat implicit-solvent
chi-squared as a check that the workflow runs, never as a verdict on a force
field.

``run_md_explicit`` is the answer to that, and it needs a GPU to be practical:
solvating the same protein takes it to roughly 25-30 thousand atoms. On a
laptop CPU that is hopeless, which is why the companion Colab notebook exists.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

__all__ = ["PreparedStructure", "prepare", "sanitise_pdb", "run_md", "MDResult",
           "benchmark", "SolvatedSystem", "solvate", "run_md_explicit"]


@dataclass
class PreparedStructure:
    path: Path
    n_atoms: int
    n_residues: int
    chain: str | None
    mutations: list[str]
    sequence: str


def sanitise_pdb(pdb_path: str | Path, out_path: str | Path,
                 chain: str | None = None, altloc: str = "A") -> dict[str, int]:
    """Reduce a deposited PDB to something a force field can accept.

    Three things in real deposited files break downstream tools, all of them
    present in 2vb1 (lysozyme at 0.65 A), which is this project's reference
    structure:

    * **Existing hydrogens.** High-resolution structures resolve hydrogens, and
      re-adding them on top produces duplicates.
    * **Alternate conformations.** Partially occupied side chains appear several
      times with different altloc codes; only one set can be simulated.
    * **Repeated chain identifiers.** 2vb1 has two chains both labelled ``A``
      (protein and solvent), so selecting a chain by identifier alone keeps both.
      Filtering here, on the text, is unambiguous.

    Returns a count of what was dropped so callers can report it.
    """
    src = Path(pdb_path).read_text(errors="replace").splitlines()
    stats = {"kept": 0, "hydrogens": 0, "altloc": 0, "other_chain": 0}
    out: list[str] = []

    for line in src:
        rec = line[:6]
        if rec in ("ATOM  ", "HETATM"):
            if chain is not None and len(line) > 21 and line[21] != chain:
                stats["other_chain"] += 1
                continue
            element = line[76:78].strip().upper() if len(line) > 77 else ""
            name = line[12:16].strip()
            if element == "H" or (not element and name[:1] == "H"):
                stats["hydrogens"] += 1
                continue
            code = line[16] if len(line) > 16 else " "
            if code not in (" ", altloc):
                stats["altloc"] += 1
                continue
            out.append(line[:16] + " " + line[17:])
            stats["kept"] += 1
        elif rec.startswith(("TER", "END", "CRYST", "MODEL", "ENDMDL", "SSBOND")):
            out.append(line)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    Path(out_path).write_text("\n".join(out) + "\nEND\n")
    return stats


def prepare(pdb_path: str | Path, out_path: str | Path, chain: str | None = None,
            mutations: list[str] | None = None, ph: float = 7.0,
            keep_heterogens: bool = False) -> PreparedStructure:
    """Clean a PDB for simulation: one chain, no solvent, complete, protonated.

    ``mutations`` uses PDBFixer's notation, e.g. ``["ALA-57-GLY"]`` meaning
    residue 57 becomes glycine. PDBFixer rebuilds the side chain in a default
    rotamer without repacking neighbours, which is adequate for conservative
    substitutions after minimisation and poor for small-to-large ones.
    """
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer

    out_path = Path(out_path)
    scratch = out_path.with_name(out_path.stem + "_sanitised.pdb")
    sanitise_pdb(pdb_path, scratch, chain=chain)

    fixer = PDBFixer(filename=str(scratch))

    if mutations:
        if chain is None:
            raise ValueError("a chain must be named when applying mutations")
        fixer.applyMutations(list(mutations), chain)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=keep_heterogens)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)

    residues = [r for r in fixer.topology.residues()]
    seq = "".join(_three_to_one(r.name) for r in residues)
    return PreparedStructure(
        path=out_path,
        n_atoms=fixer.topology.getNumAtoms(),
        n_residues=len(residues),
        chain=chain,
        mutations=list(mutations or []),
        sequence=seq,
    )


@dataclass
class MDResult:
    trajectory: Path
    topology: Path
    n_frames: int
    simulated_ps: float
    wall_seconds: float
    platform: str
    ns_per_day: float

    def __repr__(self) -> str:
        return (f"MDResult({self.n_frames} frames, {self.simulated_ps:.0f} ps, "
                f"{self.platform}, {self.ns_per_day:.2f} ns/day)")


def run_md(prepared: PreparedStructure, out_dir: str | Path,
           production_ps: float = 50.0, equilibration_ps: float = 5.0,
           frame_interval_ps: float = 1.0, temperature_K: float = 300.0,
           timestep_fs: float = 4.0, implicit: str = "implicit/gbn2.xml",
           cutoff_nm: float = 1.0, hydrogen_mass_amu: float = 4.0,
           platform: str | None = None, threads: int | None = None,
           seed: int = 0) -> MDResult:
    """Minimise, equilibrate, then run production MD in implicit solvent.

    The defaults were chosen by measurement on a laptop CPU, where they give
    about an elevenfold speedup over the obvious settings:

    ``cutoff_nm=1.0``
        The dominant cost. A 2 nm cutoff encloses every pair of a protein this
        size, so it buys nothing while leaving the calculation fully O(N^2).
        Shortening it is worth roughly 4x. Measured radius of gyration is
        unchanged across 2.0, 1.5, 1.2 and 1.0 nm (14.21-14.25 A, within the
        thermal fluctuation), so for a compact globular protein the saving is
        free. For an extended or multidomain system, check before trusting it.

    ``timestep_fs=4.0`` with ``hydrogen_mass_amu=4.0``
        Hydrogen mass repartitioning moves mass onto hydrogens so the fastest
        bond vibrations no longer set the stable timestep. Worth about 2x, and
        standard practice. The two belong together: a 4 fs step without the
        repartitioning is unstable.

    ``threads=None``
        Leave OpenMM to choose. Pinning the count measured slower in every
        combination tried.
    """
    import openmm as mm
    from openmm import app, unit

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb = app.PDBFile(str(prepared.path))
    forcefield = app.ForceField("amber14-all.xml", implicit)
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=app.HBonds,
        hydrogenMass=hydrogen_mass_amu * unit.amu,
        soluteDielectric=1.0,
        solventDielectric=78.5,
    )

    integrator = mm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin, 1.0 / unit.picosecond,
        timestep_fs * unit.femtosecond)
    integrator.setRandomNumberSeed(seed)

    plat = mm.Platform.getPlatformByName(platform) if platform else None
    props = {"Threads": str(threads)} if (plat and threads) else None
    if plat is not None:
        sim = app.Simulation(pdb.topology, system, integrator, plat, props)
    else:
        sim = app.Simulation(pdb.topology, system, integrator)
    sim.context.setPositions(pdb.positions)

    sim.minimizeEnergy()
    sim.context.setVelocitiesToTemperature(temperature_K * unit.kelvin, seed)

    steps_per_ps = 1000.0 / timestep_fs
    sim.step(int(equilibration_ps * steps_per_ps))

    traj_path = out_dir / "production.dcd"
    interval = max(int(frame_interval_ps * steps_per_ps), 1)
    sim.reporters.append(app.DCDReporter(str(traj_path), interval))

    n_steps = int(production_ps * steps_per_ps)
    t0 = time.perf_counter()
    sim.step(n_steps)
    wall = time.perf_counter() - t0

    simulated_ns = production_ps / 1000.0
    return MDResult(
        trajectory=traj_path,
        topology=prepared.path,
        n_frames=n_steps // interval,
        simulated_ps=production_ps,
        wall_seconds=wall,
        platform=sim.context.getPlatform().getName(),
        ns_per_day=simulated_ns / (wall / 86400.0) if wall > 0 else float("nan"),
    )


def benchmark(prepared: PreparedStructure, ps: float = 10.0,
              platforms: list[str] | None = None) -> dict[str, float]:
    """Measure ns/day for each available platform.

    Worth running once before planning anything: throughput on an integrated
    GPU or a laptop CPU varies by more than an order of magnitude, and every
    scheduling decision downstream depends on the real number.
    """
    import openmm as mm
    import tempfile

    if platforms is None:
        platforms = [mm.Platform.getPlatform(i).getName()
                     for i in range(mm.Platform.getNumPlatforms())]

    results: dict[str, float] = {}
    for name in platforms:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                res = run_md(prepared, tmp, production_ps=ps, equilibration_ps=2.0,
                             frame_interval_ps=ps, platform=name)
            results[name] = res.ns_per_day
        except Exception as exc:
            results[name] = float("nan")
            print(f"  {name}: unavailable ({type(exc).__name__}: {exc})")
    return results


@dataclass
class SolvatedSystem:
    """A protein in an explicit water box, ready to simulate."""

    path: Path                 # full solvated system
    protein_pdb: Path          # protein only, the topology for analysis
    n_atoms: int
    n_waters: int
    n_protein_atoms: int
    box_nm: tuple[float, float, float]

    def __repr__(self) -> str:
        return (f"SolvatedSystem({self.n_atoms} atoms, {self.n_waters} waters, "
                f"box {self.box_nm[0]:.1f}x{self.box_nm[1]:.1f}x{self.box_nm[2]:.1f} nm)")


def _protein_indices(topology) -> list[int]:
    """Atom indices that are neither solvent nor counter-ions."""
    solvent = {"HOH", "WAT", "NA", "CL", "K", "MG", "CA", "ZN", "SOD", "CLA"}
    return [a.index for a in topology.atoms() if a.residue.name.upper() not in solvent]


def _describe_solvated(path: Path) -> SolvatedSystem:
    """Rebuild the handle for an already-written solvated system."""
    from openmm import app, unit

    pdb = app.PDBFile(str(path))
    box = pdb.topology.getUnitCellDimensions().value_in_unit(unit.nanometer)
    return SolvatedSystem(
        path=path,
        protein_pdb=path.with_name(path.stem + "_protein.pdb"),
        n_atoms=pdb.topology.getNumAtoms(),
        n_waters=sum(1 for r in pdb.topology.residues()
                     if r.name.upper() in ("HOH", "WAT")),
        n_protein_atoms=len(_protein_indices(pdb.topology)),
        box_nm=(box[0], box[1], box[2]))


def solvate(prepared: PreparedStructure, out_path: str | Path,
            padding_nm: float = 1.0, ionic_strength_M: float = 0.15,
            water_model: str = "tip3p", reuse: bool = True,
            forcefield_files: tuple[str, ...] = ("amber14-all.xml",)) -> SolvatedSystem:
    """Put the prepared protein in a neutralised, salted water box.

    ``padding_nm`` is the minimum distance from the protein to the box edge.
    One nanometre is the usual choice: enough that the protein does not
    interact with its own periodic image under a 1 nm cutoff, without paying
    for water that contributes nothing.

    A protein-only PDB is written alongside the solvated system. Trajectories
    from ``run_md_explicit`` store protein coordinates only -- water dominates
    the atom count but contributes nothing the SAXS calculator uses, and
    keeping it would multiply file sizes by an order of magnitude -- so that
    file is the topology the analysis needs.

    With ``reuse`` set, an existing output is described and returned rather
    than rebuilt. Adding water is single-threaded and slow -- 48 minutes for
    this protein on the laptop used to develop the project -- so repeating it
    after a disconnected cloud session is worth avoiding. Pass ``reuse=False``
    when changing padding or ionic strength, since the cached file would
    otherwise silently win.
    """
    from openmm import app, unit

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if reuse and out_path.exists() and out_path.with_name(
            out_path.stem + "_protein.pdb").exists():
        return _describe_solvated(out_path)

    pdb = app.PDBFile(str(prepared.path))
    water_xml = f"amber14/{water_model}.xml"
    forcefield = app.ForceField(*forcefield_files, water_xml)

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(forcefield,
                        model=water_model,
                        padding=padding_nm * unit.nanometer,
                        ionicStrength=ionic_strength_M * unit.molar,
                        neutralize=True)

    with open(out_path, "w") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)

    keep = _protein_indices(modeller.topology)
    protein_pdb = out_path.with_name(out_path.stem + "_protein.pdb")
    sub = app.Modeller(modeller.topology, modeller.positions)
    sub.delete([a for a in sub.topology.atoms() if a.index not in set(keep)])
    with open(protein_pdb, "w") as fh:
        app.PDBFile.writeFile(sub.topology, sub.positions, fh, keepIds=True)

    box = modeller.topology.getUnitCellDimensions().value_in_unit(unit.nanometer)
    waters = sum(1 for r in modeller.topology.residues()
                 if r.name.upper() in ("HOH", "WAT"))
    return SolvatedSystem(path=out_path, protein_pdb=protein_pdb,
                          n_atoms=modeller.topology.getNumAtoms(),
                          n_waters=waters, n_protein_atoms=len(keep),
                          box_nm=(box[0], box[1], box[2]))


def run_md_explicit(solvated: SolvatedSystem, out_dir: str | Path,
                    production_ns: float = 2.0, equilibration_ps: float = 200.0,
                    frame_interval_ps: float = 10.0, temperature_K: float = 300.0,
                    pressure_bar: float = 1.0, timestep_fs: float = 4.0,
                    hydrogen_mass_amu: float = 4.0, cutoff_nm: float = 1.0,
                    water_model: str = "tip3p", platform: str | None = None,
                    seed: int = 0, progress: bool = True,
                    checkpoint: str | Path | None = None) -> MDResult:
    """Explicit-solvent MD under PME, NPT equilibration then NPT production.

    Practical only on a GPU. A solvated small protein is roughly 25-30 thousand
    atoms, where a laptop CPU manages well under a nanosecond per day; a modest
    GPU does one to two orders of magnitude better.

    The trajectory stores protein atoms only. Water is most of the system but
    none of it reaches the SAXS calculator, which models the hydration shell
    implicitly through its fitted contrast parameter.

    ``checkpoint`` makes the run survive interruption, which matters on a free
    cloud runtime that can disconnect mid-simulation. If the file exists the
    run resumes from it and skips minimisation and equilibration; otherwise it
    starts fresh and writes to it periodically. Put it somewhere that outlives
    the session -- on Colab that means Google Drive, not the local disk.

    Frames from before an interruption are kept. Each call writes its own
    ``segment###.dcd`` and they are joined into one trajectory at the end,
    because DCD cannot be appended to -- see ``_merge_segments``. Resuming is
    idempotent in the sense that asking for a length already reached does
    nothing and leaves the existing frames alone.
    """
    import mdtraj as mdt
    import openmm as mm
    from mdtraj.reporters import DCDReporter as SubsetDCDReporter
    from openmm import app, unit

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb = app.PDBFile(str(solvated.path))
    forcefield = app.ForceField("amber14-all.xml", f"amber14/{water_model}.xml")
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=app.HBonds,
        hydrogenMass=hydrogen_mass_amu * unit.amu,
        rigidWater=True,
    )
    system.addForce(mm.MonteCarloBarostat(pressure_bar * unit.bar,
                                          temperature_K * unit.kelvin))

    integrator = mm.LangevinMiddleIntegrator(
        temperature_K * unit.kelvin, 1.0 / unit.picosecond,
        timestep_fs * unit.femtosecond)
    integrator.setRandomNumberSeed(seed)

    plat = mm.Platform.getPlatformByName(platform) if platform else None
    sim = (app.Simulation(pdb.topology, system, integrator, plat) if plat
           else app.Simulation(pdb.topology, system, integrator))
    sim.context.setPositions(pdb.positions)

    steps_per_ps = 1000.0 / timestep_fs
    interval = max(int(frame_interval_ps * steps_per_ps), 1)
    n_steps = int(production_ns * 1000.0 * steps_per_ps)
    traj_path = out_dir / "production_protein.dcd"
    keep = _protein_indices(pdb.topology)

    checkpoint = Path(checkpoint) if checkpoint else None
    resumed = checkpoint is not None and checkpoint.exists()
    if resumed:
        sim.loadCheckpoint(str(checkpoint))
        done = sim.context.getStepCount()
        n_steps = max(n_steps - done, 0)
        print(f"resumed from {checkpoint} at production step {done} "
              f"({done * timestep_fs / 1000.0:.1f} ps); "
              f"{n_steps} steps remaining", flush=True)
    else:
        done = 0
        sim.minimizeEnergy()
        sim.context.setVelocitiesToTemperature(temperature_K * unit.kelvin, seed)
        sim.step(int(equilibration_ps * steps_per_ps))
        # Zero the counter so the checkpoint's step count means "production
        # steps completed". Without this the equilibration steps are counted
        # against the production budget on resume, and a resumed run stops
        # early -- or, if equilibration alone exceeds the budget, does nothing.
        sim.context.setStepCount(0)

    # Each run writes its own segment, merged into one trajectory at the end.
    # The DCD format cannot be appended to -- mdtraj opens it 'r' or 'w' only --
    # so writing straight to the final path would discard every frame recorded
    # before an interruption, which defeats the point of checkpointing.
    if not resumed:
        # A fresh start is a new trajectory. Leaving segments from an earlier
        # run in place would silently splice two unrelated simulations
        # together at merge time.
        for stale in out_dir.glob("segment*.dcd"):
            stale.unlink()
    segment = out_dir / f"segment{len(list(out_dir.glob('segment*.dcd'))):03d}.dcd"
    sim.reporters.append(SubsetDCDReporter(str(segment), interval, atomSubset=keep))
    if checkpoint is not None:
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        sim.reporters.append(app.CheckpointReporter(str(checkpoint), interval))
    if progress:
        # totalSteps counts from the start of production, not from this call,
        # so the remaining-time estimate stays right across a resume.
        sim.reporters.append(app.StateDataReporter(
            __import__("sys").stdout, max(interval * 5, 1), step=True,
            potentialEnergy=True, temperature=True, speed=True,
            remainingTime=True, totalSteps=done + n_steps))

    t0 = time.perf_counter()
    sim.step(n_steps)
    wall = time.perf_counter() - t0

    for reporter in list(sim.reporters):
        close = getattr(reporter, "close", None)
        if close is not None:
            try:
                close()
            except Exception:
                pass
    sim.reporters.clear()

    n_frames = _merge_segments(out_dir, solvated.protein_pdb, traj_path)

    ns_done = n_steps * timestep_fs / 1e6          # this call only, for the rate
    return MDResult(
        trajectory=traj_path,
        topology=solvated.protein_pdb,
        n_frames=n_frames,
        # Cumulative across resumes, so a run cut short does not claim the
        # full requested length.
        simulated_ps=(done + n_steps) * timestep_fs / 1000.0,
        wall_seconds=wall,
        platform=sim.context.getPlatform().getName(),
        ns_per_day=ns_done / (wall / 86400.0) if wall > 0 else float("nan"),
    )


def _merge_segments(out_dir: Path, topology_pdb: Path, out_path: Path) -> int:
    """Join the per-run segments into one trajectory and return its frame count.

    A segment cut short by a hard kill may be unreadable; those are skipped
    with a warning rather than losing the whole run.
    """
    import mdtraj as mdt

    segments = sorted(out_dir.glob("segment*.dcd"))
    parts = []
    for seg in segments:
        try:
            parts.append(mdt.load(str(seg), top=str(topology_pdb)))
        except Exception as exc:
            print(f"  skipping unreadable segment {seg.name}: "
                  f"{type(exc).__name__}", flush=True)
    if not parts:
        raise RuntimeError(f"no readable trajectory segments in {out_dir}")

    joined = parts[0]
    for part in parts[1:]:
        joined = joined.join(part)
    joined.save_dcd(str(out_path))
    if len(parts) > 1:
        print(f"  merged {len(parts)} segments -> {joined.n_frames} frames", flush=True)
    return int(joined.n_frames)


_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C",
}


def _three_to_one(name: str) -> str:
    return _AA3_TO_1.get(name.upper(), "X")
