# Vertex Model

A 2D vertex model of epithelial/confluent tissue mechanics, written in Fortran, with a MATLAB toolkit for analysis and visualization. Cells are represented as polygons sharing vertices and edges with their neighbors; the simulation integrates vertex positions under mechanical, active, and biochemical forces, with topological remodeling (T1/T2 transitions) and cell division.

## Features

- **Core mechanics**: vertex-model energy (area + perimeter elasticity), T1 (neighbor exchange) and T2 (cell extrusion) transitions.
- **Cell proliferation**: area-threshold-triggered division, splitting a cell along its principal axis.
- **Motility**:
  - Apolar random motility — uniform, or spatially varying as a gradient or localized hotspots.
  - Gradient motility can be **Lagrangian** (assigned once, carried per vertex — the default) or **Eulerian** (continuously re-derived from each vertex's current position, via `if_motility_Eulerian`).
  - Active Brownian polar motility (ABP) and a separate simple polar-motility force.
  - Optional temporal decay of motility strength.
- **Mechanochemistry**: Rho/ROCK/Myosin reaction-diffusion dynamics coupled to active contractility (Euler or RK4 integration).
- **Boundary conditions**: free, fixed (all/top/bottom), sudden or oscillatory shear, tissue perturbation/squeeze, a "dirac comb" force protocol.
- **Restart support**: continue a run from a mid-simulation snapshot (`nrun=2`) without overwriting the original run's output.

## Repository layout

| Path | Description |
|---|---|
| `vertexmain.f90` | Main program — reads input, runs the timestep loop. |
| `allocation.f90` | Parameter/array declarations, `read_input`, `allocate_arrays`, `read_data`, `write_output`. |
| `Force.f90` | Mechanical + active (motility, ABP, polar) force calculation. |
| `Geometry.f90` | Polygon geometry (area, perimeter, boundary detection), division-splitting geometry. |
| `T1_swap.f90` / `T2_swap.f90` | Topological neighbor-exchange and cell-extrusion transitions. |
| `Proliferation.f90` | Cell division. |
| `Stress.f90` | Tissue stress tensor and shear protocols. |
| `System_Info.f90` | Energy calculation, cell centroids. |
| `array_info.f90` | Small array utilities (sort, unique). |
| `Generate_Initial_Mesh.f90` | Standalone Fortran initial-mesh generator (jittered hexagonal lattice + Voronoi tessellation) — a dependency-free alternative to `Main.m`. |
| `Main.m` | MATLAB initial-mesh generator (original method; still works, kept as reference/fallback). |
| `compile.sh` | Builds `vertexmain.exe` (and optionally `generate_initial_mesh.exe`). |
| `para_Simulation.dat` | Main simulation parameters (physics, feature flags, timestepping). |
| `para_MeshDims.dat` | Mesh/array dimensioning (`Lx`, `Ly`, array capacities). |
| `para_MeshGen.dat` | Parameters for `Generate_Initial_Mesh.f90`. |
| `v_in.dat`, `inn_in.dat`, `num_in.dat` | Initial mesh state (vertex positions, cell-vertex connectivity, vertices-per-cell). |
| `data/` | Simulation output (see below). |
| `Analysis_Codes_matlab/` | MATLAB analysis and plotting toolkit (see below). |
| `TCSH_sripts/` | tcsh batch scripts for cluster parameter sweeps and batched analysis runs. |
| `log.txt` | Running log of bug fixes, optimizations, and their verification — check here before assuming code behavior. |

## Requirements

- `gfortran` (GCC's Fortran compiler).
- MATLAB (or Octave, untested) for the analysis/plotting toolkit — not required to run the simulator itself.

## Building

```sh
sh compile.sh
```

Produces `vertexmain.exe`. `Generate_Initial_Mesh.f90`'s build line is commented out in `compile.sh` by default (see the comment there — only meant to be compiled/run for `nrun=1`, since it rewrites `para_MeshDims.dat`); uncomment it, or compile it directly, when you need a fresh mesh.

## Running a simulation

1. **Generate an initial mesh** (skip this if `v_in.dat`/`inn_in.dat`/`num_in.dat` already exist and you're restarting or reusing a mesh):
   - Fortran: edit `para_MeshGen.dat`, then build and run `generate_initial_mesh.exe`.
   - MATLAB: run `Main.m` (uses `voronoin`; slower, kept as the original/reference method).
2. **Set parameters** in `para_Simulation.dat` — physics constants, feature flags (motility, division, shear, ...), `totT`/`dt`/`it_dump`. Note this file is read **positionally** (plain sequential `read()` in `allocation.f90::read_input`) — if you add a parameter, insert the corresponding line at the exact matching position, and update `Analysis_Codes_matlab/ReadPara1Params.m`'s `names` list to match.
3. **Run**:
   ```sh
   ./vertexmain.exe
   ```
4. **Restart** a run from a mid-simulation snapshot: set `nrun=2` and `nrun2_initialTime` (the snapshot `it` to resume from) in `para_Simulation.dat`. Restart output is written under an `nrun2_` filename prefix so it never overwrites the original run's data.

## Output (`data/`)

Per-snapshot dumps (every `it_dump` steps): `v_<it>.dat`, `inn_<it>.dat`, `num_<it>.dat` (mesh state), `force_<it>.dat`, `Myosin_<it>.dat`, `cell_identity_<it>.dat`.
Whole-run summary series (rewritten periodically, safe to read while a run is in progress): `Energy.dat`, `ShearStress.dat`, `T1_count.dat`, `T2_count.dat`. `motility_store.dat` (written once, at `it=1`) records the initial per-vertex motility field. All of these get an `nrun2_` prefix for a restart run.

## Analysis toolkit (`Analysis_Codes_matlab/`)

- `RunPlotAnalysis.m` — top-level driver: set `nrun`/`itStart`/`itEnd`/`itInterval` and which panels to enable, calls `PlotAnalysis.m` for a multi-panel figure (Energy, Pressure, Force, Circularity, Q(t), MSD, cumulative T1/T2 counts, ...).
- `Movie_Code.m` (+ `MovieCode_halflatt.m`, `Movie_Code_WithMyosin.m` variants) — renders an AVI movie of the tissue, colored by any of `Force`/`Motility`/`Myosin`/`Rho`/`ROCK`/`Area`/`Perimeter`/`ShapeFactor`/`NumVertices` (see `ComputeCellColorData.m`).
- `LoadData.m` / `LoadGlobalTimeSeries.m` — read a snapshot's / the whole-run summary data from `data/`.
- `TisuePlot.m` — draws one tissue snapshot given loaded mesh + color data.
- `ReadPara1Params.m` — parses `para_Simulation.dat` into a named struct (`p.dt`, `p.totT`, `p.if_motility`, ...) — use this instead of hardcoding row numbers.
- `compute_*.m` — per-quantity computation helpers (circularity, MSD, Q(t), stress tensor, force stats) used by the panels above.

## Development

Ongoing bug fixes and optimizations are tracked in `log.txt` (file:line references, root cause, fix, and how each was verified) — check there before assuming existing behavior is correct or re-diagnosing a known issue. `Claude_Testing` is the working/review branch; `master` is the main branch.
