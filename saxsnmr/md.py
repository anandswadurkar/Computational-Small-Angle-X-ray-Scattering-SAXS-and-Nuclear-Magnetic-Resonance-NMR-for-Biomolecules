"""Structure preparation and short implicit-solvent molecular dynamics.

Implicit solvent is what makes this project run on a laptop: dropping explicit
water takes a small protein from tens of thousands of atoms to a couple of
thousand, which is the difference between nanoseconds per day and nanoseconds
per hour.

The cost is a known bias. GB-family solvent models tend to over-stabilise salt
bridges and favour over-compact structures, and the radius of gyration is
exactly what SAXS measures most sensitively. Treat chi-squared from these runs
as a check that the workflow is sound, not as a verdict on a force field. The
README says the same thing; it is worth repeating here because this is the
module where the approximation enters.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path

__all__ = ["PreparedStructure", "prepare", "sanitise_pdb", "run_md", "MDResult", "benchmark"]


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


_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C",
}


def _three_to_one(name: str) -> str:
    return _AA3_TO_1.get(name.upper(), "X")
