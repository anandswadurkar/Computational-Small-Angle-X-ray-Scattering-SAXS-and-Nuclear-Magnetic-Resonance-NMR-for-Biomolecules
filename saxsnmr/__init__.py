"""Validate short implicit-solvent MD against experimental SAXS and NMR data.

Typical use::

    from saxsnmr import Fetcher, prepare, run_md, ensemble_fit

    f = Fetcher(cache_dir="cache")
    curve = f.sasbdb_curve("SASDMJ2").trim(q_max=0.5)
    struct, provenance = f.structure_for("P00698")
    prep = prepare(struct, "work/prepared.pdb", chain="A")
    md = run_md(prep, "work", production_ps=200)
    fit = ensemble_fit(md, curve)

See the README for what the resulting chi-squared does and does not establish.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .fetch import Construct, Fetcher, SaxsCurve, StructureHit, summarise_constructs
from .md import (MDResult, PreparedStructure, benchmark, prepare, run_md,
                 sanitise_pdb)
from .nmr import ShiftComparison, best_offset, compare_to_dssp, csi_assignment, secondary_shifts
from .saxs import SaxsFit, ScatteringTables, fit, profile

__version__ = "0.1.0"

__all__ = [
    "Fetcher", "SaxsCurve", "Construct", "StructureHit", "summarise_constructs",
    "prepare", "sanitise_pdb", "run_md", "benchmark", "PreparedStructure", "MDResult",
    "ScatteringTables", "profile", "fit", "SaxsFit",
    "secondary_shifts", "csi_assignment", "compare_to_dssp", "best_offset",
    "ShiftComparison",
    "frame_tables", "ensemble_fit", "radius_of_gyration", "load_trajectory",
]

# Shrake-Rupley gives absolute SASA; dividing by a per-element maximum turns it
# into the fractional accessibility the hydration term expects.
_MAX_SASA = {"H": 0.6, "C": 1.2, "N": 1.1, "O": 1.1, "S": 1.5, "P": 1.5}


def load_trajectory(path: str | Path, top: str | Path | None = None):
    """Load a PDB or DCD into an ``mdtraj.Trajectory``.

    A convenience wrapper that pairs a DCD with its topology, so callers do not
    have to remember that a DCD carries no atom names and converts coordinates
    from angstrom. ``mdtraj.load`` would serve equally well; this simply keeps
    the call sites uniform across the two formats the project writes.
    """
    from mdtraj.core.trajectory import Trajectory
    from mdtraj.formats import DCDTrajectoryFile, PDBTrajectoryFile

    path = str(path)
    if path.lower().endswith(".pdb"):
        with PDBTrajectoryFile(path) as fh:
            return Trajectory(xyz=np.asarray(fh.positions) / 10.0, topology=fh.topology)

    if top is None:
        raise ValueError("a topology (.pdb) is required to load a DCD trajectory")
    with PDBTrajectoryFile(str(top)) as fh:
        topology = fh.topology
    with DCDTrajectoryFile(path) as fh:
        xyz, _lengths, _angles = fh.read()
    return Trajectory(xyz=np.asarray(xyz) / 10.0, topology=topology)


def _elements_and_sasa(traj, frame: int):
    import mdtraj as md

    sasa = md.shrake_rupley(traj[frame], mode="atom")[0] * 100.0  # nm^2 -> A^2
    elements, frac = [], []
    for atom, s in zip(traj.topology.atoms, sasa):
        sym = atom.element.symbol if atom.element is not None else "C"
        if sym not in _MAX_SASA:
            sym = "C"
        elements.append(sym)
        frac.append(min(s / (4.0 * np.pi * _MAX_SASA[sym] ** 2), 1.0))
    return elements, np.asarray(frac)


def frame_tables(traj, frame: int, bin_width: float = 0.25,
                 r_max: float | None = None) -> ScatteringTables:
    """Scattering tables for one trajectory frame (coordinates in angstrom)."""
    elements, frac = _elements_and_sasa(traj, frame)
    xyz = traj.xyz[frame] * 10.0
    return ScatteringTables.from_coordinates(xyz, elements, frac,
                                             bin_width=bin_width, r_max=r_max)


def radius_of_gyration(traj) -> np.ndarray:
    """Rg per frame in angstrom -- the quantity SAXS constrains most tightly."""
    import mdtraj as md

    return md.compute_rg(traj) * 10.0


def ensemble_fit(md_result: "MDResult | object", curve: SaxsCurve,
                 stride: int = 1, max_frames: int | None = None,
                 bin_width: float = 0.25) -> SaxsFit:
    """Average I(q) over trajectory frames, then fit scale, c1 and c2 once.

    Averaging intensities before fitting (rather than fitting each frame) is the
    correct order: the experiment measures one ensemble-averaged curve, so the
    model should be the ensemble average too.
    """
    traj_path = getattr(md_result, "trajectory", md_result)
    top_path = getattr(md_result, "topology", None)
    traj = load_trajectory(traj_path, top=top_path)

    frames = list(range(0, traj.n_frames, stride))
    if max_frames is not None:
        frames = frames[:max_frames]
    if not frames:
        raise ValueError("trajectory contained no frames to average")

    # One shared r_max so every frame bins identically and the tables can be
    # averaged. Taking it from the whole trajectory's extent avoids a frame
    # that happens to expand slightly falling outside the range.
    xyz = traj.xyz[frames] * 10.0
    span = float(np.max([np.ptp(f, axis=0).max() for f in xyz]))
    r_max = span * 1.75 + 10.0

    tables = [frame_tables(traj, i, bin_width=bin_width, r_max=r_max) for i in frames]
    mean_tables = ScatteringTables.average(tables)
    return fit(mean_tables, curve.q, curve.intensity, curve.sigma)
