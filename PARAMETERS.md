# Parameter reference

`para_MeshGen.dat` and `para_Simulation.dat` are both read **positionally** (plain
sequential `read()` calls — `Generate_Initial_Mesh.f90` and
`allocation.f90::read_input` respectively), so their inline comments are kept short
— just enough to identify each line, not to explain it. This file documents every
parameter in both, in the exact order they're read, with a little more explanation
where the one-line comment isn't the whole story. If a parameter's inline comment
already says everything there is to say, its entry here is just as short.

If you add a parameter complex enough to need more than a few words, give it the
same treatment: a short pointer in the `.dat` file (`-- see PARAMETERS.md`), the
real explanation here — and remember both files must stay positionally in sync
with their readers (`Generate_Initial_Mesh.f90`, `allocation.f90::read_input`,
`Analysis_Codes_matlab/ReadPara1Params.m`).

## `para_MeshGen.dat` — read by `Generate_Initial_Mesh.f90`

| Parameter | Description |
|---|---|
| `Lx` | Lattice columns. Also doubles as the exact physical box width once the mesh is built (lattice spacing is 1) — used that way everywhere else in the codebase (`mot_Lc*Ly`, the stress-tensor core radius, `Lx_box`/`Ly_box` under PBC). |
| `Ly` | Lattice rows. Must be even (the hexagonal-offset lattice alternates row offset by parity). Same dual role as `Lx`. |
| `seed` | RNG seed for the lattice jitter. Fixed value → the exact same mesh every time; change it for a different random mesh. |
| `jitter_interior_x` | Max random x-offset applied to interior-cell seed points, as a fraction of lattice spacing. |
| `jitter_interior_y` | Same, y-offset. One random draw per point scales *both* its x and y jitter (correlated direction, not independent noise) — see the module header comment. |
| `jitter_boundary_x` | Max random x-offset for outer-ring seed points specifically. Normally smaller than the interior jitter to keep a free mesh's edge less ragged. Ignored (every cell uses the interior jitter instead) when `if_periodic=.true.` — see that entry. |
| `jitter_boundary_y` | Same, y-offset for the outer ring. |
| `cell_headroom` | Extra cell slots reserved beyond `Lx*Ly`, for cells created later by division. Sets `num_dim = Lx*Ly + cell_headroom` in the written `para_MeshDims.dat`. |
| `vertex_slot_headroom` | Extra per-cell vertex-list slots reserved beyond the largest polygon actually generated. Sets `inn_dim1`. |
| `vertex_pool_headroom` | Extra vertex-array rows reserved beyond the vertex count actually generated (new vertices come from T1 flips and division). Sets `v_dim2`; rule of thumb ~2×N. |
| `if_periodic` | `.false.` (default) — the ordinary free-boundary mesh every earlier feature assumes: ghost seed copies are used only to clip boundary-cell shapes against something, not to connect them to anything (see the module header — "a boundary-conditioning trick, not true periodic BCs"). `.true.` — a genuinely wrap-connected torus mesh: a post-pass identifies vertices independently computed by cells on opposite edges of the box (their raw positions differ by exactly one lattice period) and unifies their IDs via union-find, so opposite edges share real topology. Verified via Euler characteristic (V−E+F=0, the torus signature), uniform vertex degree 3, and exact area conservation (`log.txt`). Required (and checked at simulation startup, `allocation.f90::read_data`) if you're going to run with `if_PBC=.true.` in `para_Simulation.dat`. |

## `para_Simulation.dat` — read by `allocation.f90::read_input`

### General / restart

| Parameter | Description |
|---|---|
| `nrun` | `1` = fresh start from `v_in.dat`/`inn_in.dat`/`num_in.dat`. `2` = restart from a mid-run snapshot at `nrun2_initialTime`; all restart output is written under an `nrun2_` filename prefix so it never overwrites the original run's data. |
| `nrun2_initialTime` | The snapshot `it` to resume from when `nrun=2`. Ignored when `nrun=1`. |
| `Ao` | Target (preferred) cell area in the area-elasticity energy term. |
| `Co` | Target (preferred) cell perimeter in the perimeter-elasticity energy term. |
| `lambda` | Area-elasticity modulus. |
| `beta` | Perimeter-elasticity modulus. Non-dimensionalized once at startup (divided by `lambda*Ao`) — the value here is the physical one, not the rescaled internal one. Silently overridden every step if `if_active_contractility` or `if_RhoROCK` is on (the two aren't mutually exclusive in code; if both are on, `if_active_contractility` wins — an unresolved design gap, `log.txt`). |
| `gamm` | Line-tension modulus (linear-in-perimeter energy term). Non-dimensionalized once at startup like `beta`. |
| `eta` | Damping/mobility coefficient — the sole friction constant relating force to velocity (`v += dt*F/eta`); uniform across the whole tissue, no per-vertex or per-cell variation exists. |
| `totT` | Total number of timesteps. |
| `dt` | Timestep size. |

### T1 / T2 topological transitions

| Parameter | Description |
|---|---|
| `if_Do_T1` | Enable T1 neighbor-exchange transitions (short-edge flips). |
| `min_d_T1` | Edge-length threshold below which an edge is a T1 candidate. |
| `if_Do_T2` | Enable T2 cell-extrusion transitions (triangular-cell removal). |
| `min_area_T2` | Area threshold below which a 3-vertex (triangular) cell is extruded. Only 3-vertex cells can T2 — no support for extruding higher-order polygons. |

### Boundary conditions

| Parameter | Description |
|---|---|
| `if_Fixed_boundary` | Pin (zero the force on) every vertex on the tissue's free outer edge. |
| `if_bottom_borders_fixed` | Pin only the bottom-most border vertices. Independent of `if_Fixed_boundary` — can be used alone for a "one wall" setup. |
| `if_top_borders_fixed` | Pin only the top-most border vertices. |
| `if_PBC` | Periodic boundary conditions — periodic in both x and y (a torus). See `Generate_Initial_Mesh.f90`'s `if_periodic` above for the matching mesh requirement, and the note there on the safety check at startup. Mutually exclusive with all three flags above (no free edge exists to pin under full periodicity), with `if_squeeze_tissue`/`if_limb_force` (both assume a real edge/corners), with `if_Shear_tissue` (today's shear mechanism directly displaces every vertex, not equivalent to genuine periodic/Lees-Edwards shear — not yet implemented), and with `if_motility_gradient`+`if_motility_Eulerian` together (that gradient formula isn't periodic in y). All checked at startup (`read_input` stops with an error if violated). Full design rationale and verification: `log.txt`. |

### Output cadence

| Parameter | Description |
|---|---|
| `it_dump` | Write a full mesh snapshot (`v`/`inn`/`num`/`force`/`Myosin`/`cell_identity`) every this many steps (always also at `it=1,2`). Whole-run summary series (`Energy`/`ShearStress`/`T1_count`/`T2_count`) are rewritten periodically at a coarser, automatically-derived cadence (`summary_dump_interval`, capped to ~200 rewrites over the whole run) — no separate parameter needed. |
| `T1_time_interval` | How often (in units of `it`) the T1 attempt loop runs. |
| `T2_time_interval` | Declared but dead — `Do_T2` is always called from inside the `T1_time_interval` loop regardless of this value (`log.txt`). |

### Motility

| Parameter | Description |
|---|---|
| `if_motility` | Enable apolar random motility (an independent random force per vertex each step, magnitude set by the per-vertex field `mot`). |
| `if_motility_gradient` | `.false.` — `mot` is uniform, equal to `etas_max` everywhere. `.true.` — `mot` varies with y (see `etas_max`/`mot_Lc`); *when* it's evaluated is controlled by `if_motility_Eulerian`. |
| `if_motility_Eulerian` | Only matters when `if_motility_gradient=.true.`. `.false.` (Lagrangian, default) — `mot` is assigned once at `it=1` from each vertex's position then, and carried by that vertex wherever it moves afterward (a persistent tag, not a standing field — measurably erodes over a long run, mainly from cell division creating new vertices whose motility must be estimated, not from passive drift). `.true.` (Eulerian) — re-derived from each vertex's *current* position every step, so the gradient stays spatially anchored regardless of tissue motion/mixing/division (~6% overhead). Not periodic in y — disallowed together with `if_PBC`. Full verification (seeded A/B): `log.txt`. |
| `etas_max` | Apolar motility strength — the uniform value when `if_motility_gradient=.false.`, or the peak value (at y=0) when the gradient is on, or the peak value at a hotspot center when `if_motility_hotspot` is on. |
| `etas_min` | Read but never actually used in the gradient formula (only `etas_max`/`mot_Lc` matter) — vestigial. |
| `mot_Lc` | Motility-gradient decay length scale, as a fraction of `Ly`: `mot = etas_max*exp(-y/(mot_Lc*Ly))`. |
| `if_motility_decay` | Exponentially decay `mot` toward 0 over `motility_decay_timeScale`. Largely moot combined with `if_motility_Eulerian` — the next step's recompute discards the decayed value before it can accumulate. |
| `motility_decay_timeScale` | Decay time constant, in units of `dt`. |
| `if_motility_hotspot` | Enable localized Gaussian motility enhancement around 1, 2, or 4 fixed hotspot cells (hardcoded cell indices near mesh corners, not user-configurable positions), instead of the y-gradient. Reuses `etas_max` as the peak value. Assigned once at `it=1` (or at restart), like Lagrangian-mode gradient — not re-evaluated every step. |
| `number_of_hotspot` | Must be 1, 2, or 4 — anything else stops with an error. |
| `sigma_hotspot` | Gaussian width of each hotspot's motility enhancement. |

### Shear

| Parameter | Description |
|---|---|
| `if_Shear_tissue` | Master switch — must be on for either shear mode below (and for `Calculate_StressTensor` to run at all) to have any effect. |
| `if_Sudden_Shearing` | One-time affine shear (`x += sudden_shearStrength*y`) applied at `it==sudden_shearWhen`. |
| `sudden_shearStrength` | Shear strain applied by the sudden-shear event. |
| `sudden_shearWhen` | Timestep the sudden shear is applied at. |
| `if_Oscillatory_Shearing` | Continuous oscillatory affine shear, applied every step once `it > Oscl_shearWhen`. |
| `Oscl_shearStrength` | Oscillatory shear strain amplitude. |
| `Oscl_shearWhen` | Timestep oscillatory shearing begins. |
| `Oscl_freq_wo` | Oscillatory shear angular frequency. |

### Tissue perturbation

| Parameter | Description |
|---|---|
| `if_perturb_tissue` | Master switch for the sinusoidal perturbation below. |
| `if_sin_perturb` | Apply a sinusoidal x-displacement, `dt*sin_perturb_strength*sin(2*pi*sin_perturb_waveNumber*y/Ly)`, to every vertex. **Note:** this currently does nothing at all unless `if_dirac_comb` is *also* `.true.` — the displacement is nested inside the comb gate (`Force.f90`), likely an unintended coupling between two conceptually independent flags (`log.txt`). |
| `sin_perturb_when` | Timestep the perturbation begins (unit of `it`). |
| `sin_perturb_strength` | Perturbation displacement amplitude. |
| `sin_perturb_waveNumber` | Number of spatial oscillations across the tissue height. |
| `if_dirac_comb` | Gate the perturbation on/off in a square-wave pulse train (see `if_sin_perturb`'s note — currently required for any sinusoidal perturbation to apply at all, even a continuous one, via a long `comb_onPeriod`). |
| `comb_onPeriod` | Pulse "on" duration, unit of `it`. |
| `comb_offPeriod` | Pulse "off" duration, unit of `it`. |

### Limb force

| Parameter | Description |
|---|---|
| `if_limb_force` | Apply a fixed outward x-force at four hardcoded near-corner cell locations (pushes left-half cells left, right-half cells right) — assumes a bounded domain with real corners; incompatible with `if_PBC`. |
| `limb_force_strength` | Force magnitude for the above. |

### Cell division

| Parameter | Description |
|---|---|
| `if_cell_division` | Enable area-threshold-triggered division (split along the principal axis). Daughter cells correctly inherit split-vertex motility; `Rho`/`ROCK`/`Myosin`/`cell_identity` are **not** propagated to daughters (`log.txt`) — a known gap for any lineage-tracking or RhoROCK+division study. |
| `area_0` | Area threshold above which a cell is a division candidate. |

### Squeeze

| Parameter | Description |
|---|---|
| `if_squeeze_tissue` | One-time uniform compression of every vertex's y-coordinate toward the current bottom border, at `it==squeeze_when`. Assumes a real bottom edge to anchor to — incompatible with `if_PBC`. |
| `squeeze_when` | Timestep the squeeze is applied at. |
| `percent_squeeze` | Compression amount, percent. |

### Active contractility

| Parameter | Description |
|---|---|
| `if_active_contractility` | Ornstein-Uhlenbeck fluctuating contractility: a per-cell "Myosin" level relaxing to the basal `beta` with correlation time `tau_contr`, directly setting the perimeter-elasticity modulus `beta` per cell. Simpler and independent of the full RhoROCK system below (see `beta`'s entry for the interaction if both are on). |
| `tau_contr` | Correlation time for the OU process. |
| `active_contr_strength` | Noise strength for the OU process. |

### Active Brownian motility (ABP)

| Parameter | Description |
|---|---|
| `if_ABP` | Self-propelled motion: each **cell** (not vertex) has one orientation, rotationally diffusing, with a propulsion force `(vo*cos(theta), vo*sin(theta))` applied identically to every vertex of that cell. No neighbor-alignment — orientations diffuse completely independently per cell (no collective/flocking behavior emerges from this alone). |
| `vo` | Self-propulsion speed. |
| `Dr` | Rotational diffusion constant for the orientation angle. |

### RhoROCK mechanochemical regulation

| Parameter | Description |
|---|---|
| `if_RhoROCK` | Per-cell Rho→ROCK→Myosin reaction cascade (area-dilation-driven Rho activation via a Hill function, linear production/decay down the chain), feeding back into the perimeter-elasticity modulus `beta` via `Myosin_Coupling_Strength`. **Despite the naming, `D_rho`/`D_ROCK`/`D_Myosin` are pure local decay rates, not diffusion** — there is no neighbor-coupling term anywhere, so Rho/ROCK/Myosin never actually spread between cells (`log.txt`). See `beta`'s entry for the interaction with `if_active_contractility` if both are on. |
| `if_RK4` | Integrate the RhoROCK ODEs with 4th-order Runge-Kutta instead of Euler. |
| `nhill` | Hill-function exponent for area-driven Rho activation. |
| `K_hill` | Hill-function half-saturation constant. |
| `A_Rho` | Rho production rate constant. |
| `A_ROCK` | ROCK production rate constant (from Rho). |
| `A_Myosin` | Myosin production rate constant (from ROCK). |
| `D_rho` | Rho decay rate (not diffusion — see `if_RhoROCK`). |
| `D_ROCK` | ROCK decay rate (not diffusion). |
| `D_Myosin` | Myosin decay rate (not diffusion). |
| `Myosin_Coupling_Strength` | Scales Myosin's feedback into the perimeter-elasticity modulus `beta`. |
| `if_myosin_noise` | Add noise directly to the Myosin ODE. |
| `if_gaussian_noise` | Use Gaussian (vs. uniform) noise when `if_myosin_noise` is on. |
| `Myosin_noise_strength` | Myosin noise magnitude. |
| `if_coupling_noise` | Add noise to the Myosin→`beta` coupling itself, instead of/in addition to the Myosin ODE. |
| `coupling_noise_strength` | Coupling-noise magnitude. |

### Polar motility

| Parameter | Description |
|---|---|
| `if_polar_motility` | A random 2D kick per cell per step, applied to all its vertices. Despite the name, this has **no persistence/direction memory** (redrawn independently every step, unlike `if_ABP`'s diffusing orientation) — it's cell-scale white noise, not a polarized/directed active force. |
| `polar_motility_strength` | Kick magnitude. |
