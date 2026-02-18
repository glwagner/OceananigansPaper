# Benchmark Plan for Oceananigans Paper Listings

**Purpose:** Provide wall-clock timing and throughput numbers for the examples shown
in the paper, in response to Reviewer #4 (comment #473) requesting performance
data (wall-clock duration / SYPD) for the code listings.

**Baseline hardware:** Single NVIDIA H100 80 GB SXM GPU, unless otherwise noted.

**General methodology:**

1. Each benchmark script should be derived from the corresponding paper listing
   or repository script, with minimal modifications (shorter runs, disabled I/O).
2. Exclude Julia compilation time by running a short warm-up phase (a few time
   steps) before starting the timing loop.
3. Disable all output writers during the timed portion to isolate compute cost.
4. Report timing as the median of 3 independent runs.
5. Collect GPU memory usage via `CUDA.memory_status()` at peak.

---

## Listing 1 -- "Hello, ocean" (2D turbulence)

**Paper configuration:**
- `NonhydrostaticModel` on a 256x256 `RectilinearGrid` (Periodic, Periodic, Flat)
- `WENO(order=9)`, domain `(0, 2pi)^2`
- `dt = 0.01`, `stop_time = 10` (1000 time steps)
- Architecture: GPU

**What to measure:**
- Total wall-clock time for the 1000 time steps (excluding compilation)
- Time per time step (ms)
- Peak GPU memory

**Benchmark simplifications:**
- Remove CairoMakie visualization code
- Remove all output writers
- Run on GPU (paper says CPU is also fine for this size; benchmark both)

**Metrics to report:**
| Metric | Unit |
|---|---|
| Wall clock (1000 steps) | seconds |
| Time per step | ms |
| Peak GPU memory | MB |
| Peak CPU memory (CPU run) | MB |

**Notes:** This is a very small problem. The benchmark mainly demonstrates
compilation and launch overhead rather than raw throughput. Include a CPU
comparison to show that GPU overhead dominates at this scale.

---

## Listing 2 -- Passive tracer with moving source

**Paper configuration:**
- Same as Listing 1, plus one passive tracer `c` with a `circling_source`
  forcing function
- `WENO(order=9)`, `tracers=:c`, `forcing=(; c=circling_source)`
- `dt = 0.01`, runs to `t = 2.5` then to `t = 10`

**What to measure:**
- Wall-clock time for the full run (1000 time steps to t=10)
- Overhead of adding one forced tracer vs Listing 1

**Benchmark simplifications:**
- Same as Listing 1 (no visualization, no I/O)
- Time the full 1000-step segment to t=10

**Metrics to report:**
| Metric | Unit |
|---|---|
| Wall clock (1000 steps) | seconds |
| Time per step | ms |
| Overhead vs Listing 1 | % |
| Peak GPU memory | MB |

---

## Listing 3 -- AbstractOperations / lazy expression trees

**Paper configuration:**
- Diagnostic computations: vorticity `zeta = dx(v) - dy(u)`, speed
  `s = sqrt(u^2 + v^2)`, enstrophy integral `Z = Integral(zeta^2, dims=1)`
- Builds `Field` objects from AbstractOperations and calls `compute!`

**What to measure:**
- Time to `compute!` a vorticity Field on the 256x256 grid
- Time to `compute!` the speed Field
- Time to `compute!` the enstrophy integral (reduction)

**Benchmark simplifications:**
- Use the model state from a completed Listing 1 run (or set fields to random data)
- Time each `compute!` call independently, 100 repetitions each, report median

**Metrics to report:**
| Metric | Unit |
|---|---|
| compute!(vorticity) | microseconds |
| compute!(speed) | microseconds |
| compute!(enstrophy integral) | microseconds |

**Notes:** This listing is about API expressiveness, not time stepping. The
benchmark quantifies the overhead of the AbstractOperations machinery.

---

## Listing 4 -- Flow past a cylinder (DNS and LES)

**Paper configuration (DNS, Re=1000):**
- `RectilinearGrid(GPU())`, size `(2*Ny, Ny)` = `(4096, 2048)` for Re=1000
- `ImmersedBoundaryGrid` with `GridFittedBoundary(cylinder)`
- `ScalarDiffusivity(nu=1/Re)`, `Centered(order=2)`, no-slip BCs
- `stop_time = 200`, adaptive time stepping (CFL = 1.0)
- `ConjugateGradientPoissonSolver` with FFT-based preconditioner

**Paper configuration (LES, Re -> infinity):**
- Size `(6144, 2048)`, `WENO(order=9)`, no closure
- Drag law BCs on immersed boundary
- `stop_time = 200`

**What to measure:**
- Wall-clock time for a fixed number of time steps (e.g., 500 steps)
- Time per time step (includes CG pressure solve)
- Number of CG iterations per step (for the CG solver variant)
- Peak GPU memory

**Benchmark simplifications:**
- Run for 500 time steps instead of full `stop_time = 200`
- Disable JLD2 output writers
- Benchmark only the LES configuration (Re -> infinity, 6144x2048) as the
  most demanding case
- Optionally also benchmark DNS at Re=1000 (4096x2048)

**Metrics to report:**
| Metric | Unit |
|---|---|
| Wall clock (500 steps) | seconds |
| Time per step | ms |
| Avg CG iterations per step | count |
| Peak GPU memory | GB |
| Grid points updated per second | Gpts/s |

---

## Listing 5 -- Cabbeling DNS

**Paper configuration:**
- `RectilinearGrid(GPU())`, size `(4096, 1024)`, topology `(Bounded, Flat, Bounded)`
- Domain: `x = (0, 2)`, `z = (-0.5, 0)` (meters)
- `TEOS10EquationOfState`, `ScalarDiffusivity(nu=1.15e-6, kappa=(T=1e-7, S=1e-9))`
- `RungeKutta3` timestepper
- `dt = 0.2 * dx^2 / nu`, `stop_time = 100`

**What to measure:**
- Wall-clock time for a fixed number of time steps (e.g., 1000 steps)
- Time per time step
- Peak GPU memory

**Benchmark simplifications:**
- Run for 1000 time steps instead of to `stop_time = 100`
- Disable JLD2 output
- Keep RK3 timestepper (3 tendency evaluations per step)

**Metrics to report:**
| Metric | Unit |
|---|---|
| Wall clock (1000 steps) | seconds |
| Time per step | ms |
| Peak GPU memory | GB |
| Grid points updated per second | Gpts/s |

**Notes:** RK3 requires 3 sub-steps per time step, so this is more expensive per
step than QAB2. Report both time-per-step and effective throughput accounting
for the 3 sub-steps.

---

## Listing 6 -- Eady problem LES

**Paper configuration:**
- `RectilinearGrid(GPU())`, size `(1024, 1024, 64)`, topology `(Periodic, Periodic, Bounded)`
- Domain: `4 km x 4 km x 128 m` (4 m horizontal, 2 m vertical spacing)
- `WENO(order=9)`, `BackgroundFields` for geostrophic shear and buoyancy
- `FPlane(f=1e-4)`, `BuoyancyTracer`, Ri = 1
- 30 days on a single H100 GPU

**What to measure:**
- Wall-clock time for 1 simulated day
- Total wall-clock time extrapolated to 30 days
- SYPD (simulated years per day)
- Time per time step
- Peak GPU memory

**Benchmark simplifications:**
- Run for 1 simulated day (instead of 30 days) for the timed segment
- Disable output writers
- Use same grid and physics as paper

**Metrics to report:**
| Metric | Unit |
|---|---|
| Wall clock (1 sim day) | minutes |
| Extrapolated wall clock (30 sim days) | hours |
| SYPD | simulated years per day |
| Time per step (median) | ms |
| Peak GPU memory | GB |
| Grid point updates per second | Gpts/s |

**Notes:** The paper states this runs for 30 days on a single H100. Report the
actual measured SYPD for comparison with the paper's implicit claim.

---

## Listing 7 -- Tidal flow past headland (Three Tree Point)

**Paper configuration:**
- `RectilinearGrid(GPU())`, size `(384, 128, 64)` (6Nz, 2Nz, Nz with Nz=64)
- `ImmersedBoundaryGrid` with `GridFittedBottom(wedge)`
- Domain: `(-3L, 3L) x (-L, L) x (-H, 0)` with H=256m, L=1024m
- `WENO(order=9)`, `TEOS10EquationOfState`, `FPlane(latitude=47.5)`
- `PerturbationAdvectionOpenBoundaryCondition` for tidal forcing
- `stop_time = 3 days`, adaptive dt starting at 5s

**What to measure:**
- Wall-clock time for 1 simulated day
- SYPD
- Time per time step
- Peak GPU memory

**Benchmark simplifications:**
- Run for 1 simulated day instead of 3 days
- Disable all 3 JLD2 output writers (xy, xz, xyz slices)
- Disable Oceanostics `ErtelPotentialVorticity` computation (not needed for timing)

**Metrics to report:**
| Metric | Unit |
|---|---|
| Wall clock (1 sim day) | minutes |
| SYPD | simulated years per day |
| Time per step (median) | ms |
| Peak GPU memory | GB |
| Grid point updates per second | Gpts/s |

**Notes:** This is a 3D LES with open boundary conditions and immersed boundaries.
The CG Poisson solver is likely needed here. Report average CG iteration count.

---

## Listing 8 -- Tidally-forced stratified flow over seamounts

**Paper configuration:**
- `HydrostaticFreeSurfaceModel`
- `RectilinearGrid`, size `(2000, 200)`, topology `(Periodic, Flat, Bounded)`
- Domain: `x = (-1000 km, 1000 km)`, `z = (-2 km, 0)`
- `ImmersedBoundaryGrid` with `GridFittedBottom` (42 random Gaussian seamounts)
- `WENO()` for both momentum and tracer advection, `BuoyancyTracer`
- Tidal forcing at M2 frequency, `dt = 1 minute`, `stop_time = 16 * T2` (~8.3 days)

**What to measure:**
- Wall-clock time for 2 tidal periods (~24.8 hours simulated)
- SYPD
- Time per time step
- Peak GPU memory

**Benchmark simplifications:**
- Run for 2 tidal periods (~1490 minutes, ~1490 time steps at dt=1min)
  instead of full 16 periods
- Disable JLD2 output
- Use same grid and physics

**Metrics to report:**
| Metric | Unit |
|---|---|
| Wall clock (2 tidal periods) | seconds |
| SYPD | simulated years per day |
| Time per step (median) | ms |
| Peak GPU memory | GB |
| Grid point updates per second | Gpts/s |

---

## Listing 9 -- Vertical mixing parameterizations (single column)

**Paper configuration:**
- `HydrostaticFreeSurfaceModel`, single column: size 50, `z = (-200, 0)`,
  topology `(Flat, Flat, Bounded)`
- `BuoyancyTracer`, flux BCs at top for momentum and buoyancy
- Tests two closures: CATKE and k-epsilon
- `dt = 1 minute`, `stop_time = 24 hours`

**What to measure:**
- Wall-clock time for the full 24-hour simulation (1440 steps) for each closure
- Time per step for CATKE vs k-epsilon

**Benchmark simplifications:**
- None needed; this is already a very fast single-column simulation
- Run on both CPU and GPU to demonstrate overhead
- Run each closure separately

**Metrics to report:**
| Metric | Unit |
|---|---|
| Wall clock (CATKE, 24h sim) | seconds |
| Wall clock (k-epsilon, 24h sim) | seconds |
| Time per step (CATKE) | ms |
| Time per step (k-epsilon) | ms |

**Notes:** This is a 1D problem. Performance is not the point; the benchmark
demonstrates that column physics parameterizations do not add significant
overhead. Consider running this on CPU only since GPU launch overhead will
dominate for a single column.

---

## Listing 10 -- Near-global baroclinic instability on spherical grids

**Paper configuration:**
- `HydrostaticFreeSurfaceModel` on three grid types:
  `LatitudeLongitudeGrid`, `TripolarGrid`, `RotatedLatitudeLongitudeGrid`
- Size: `(4*360, 4*170, 10)` = `(1440, 680, 10)` per grid
- `WENOVectorInvariant(order=9)` momentum, `WENO(order=7)` tracers
- `TEOS10EquationOfState`, `SplitExplicitFreeSurface(substeps=60)`
- `dt = 2 minutes`, `stop_time = 180 days`

**What to measure:**
- Wall-clock time for 10 simulated days on each grid type
- SYPD for each grid type
- Peak GPU memory for each grid type

**Benchmark simplifications:**
- Run for 10 simulated days instead of 180 days
- Disable output writers
- Benchmark each grid type separately

**Metrics to report:**
| Metric | Grid | Unit |
|---|---|---|
| Wall clock (10 sim days) | LatLon / Tripolar / Rotated | minutes |
| SYPD | LatLon / Tripolar / Rotated | simulated years per day |
| Time per step (median) | LatLon / Tripolar / Rotated | ms |
| Peak GPU memory | LatLon / Tripolar / Rotated | GB |
| Grid point updates per second | LatLon / Tripolar / Rotated | Gpts/s |

**Notes:** The tripolar grid has an `ImmersedBoundaryGrid` with cylindrical
islands at the north pole singularities. This may affect performance relative
to the plain `LatitudeLongitudeGrid`. The comparison across grid types is
informative.

---

## Listing 11 -- Coupled ocean-sea ice (1/6 degree, ClimaOcean)

**Paper configuration (from sixth_degree_omip.jl and Listing 11 in PDF):**
- `TripolarGrid(GPU())`, size `(2160, 1080, 60)` -- 1/6th degree
- `ImmersedBoundaryGrid` with ETOPO1 bathymetry, `active_cells_map=true`
- `WENOVectorInvariant()` momentum, `WENO(order=7)` tracers
- `SplitExplicitFreeSurface(cfl=0.7, fixed_dt=10min)`
- CATKE vertical mixing, `CATKEVerticalDiffusivity`
- `OceanSeaIceModel` with `JRA55PrescribedAtmosphere`, sea ice model
- Spin up at `dt = 20s` for 60 days, then `dt = 6 min` for 58 years

**What to measure:**
- SYPD during the production run phase (dt = 6 minutes)
- Wall-clock time for 30 simulated days at production time step
- Peak GPU memory
- Time per time step breakdown (ocean tendency, barotropic substeps,
  sea ice, atmospheric forcing)

**Benchmark simplifications:**
- Initialize from ECCO state (requires data download)
- Run for 30 simulated days at `dt = 6 minutes` after a short 2-day spin-up
  at `dt = 20s`
- Disable all output writers and checkpointers
- Single GPU only

**Metrics to report:**
| Metric | Unit |
|---|---|
| SYPD (production dt) | simulated years per day |
| Wall clock (30 sim days) | hours |
| Time per step (median) | ms |
| Peak GPU memory | GB |
| Total ocean grid points | million |
| Active columns | million |

**Notes:** The paper (section 3.2.3, page 24) states this configuration achieves
1.5 SYPD on a single H100 GPU. The benchmark should reproduce or update this
number. This is the headline performance result for the paper. Data download
for bathymetry, ECCO, and JRA55 is required; ensure this is done before the
timed segment.

---

## Summary table of benchmark configurations

| # | Listing | Model | Grid size | Key features | Expected run time |
|---|---------|-------|-----------|-------------|-------------------|
| 1 | Hello ocean | Nonhydrostatic | 256x256 (2D) | WENO(9) | seconds |
| 2 | Circling tracer | Nonhydrostatic | 256x256 (2D) | WENO(9), forced tracer | seconds |
| 3 | AbstractOps | (diagnostic) | 256x256 (2D) | compute! benchmarks | < 1 second |
| 4 | Cylinder LES | Nonhydrostatic | 6144x2048 (2D) | ImmersedBoundary, CG solver | minutes |
| 5 | Cabbeling DNS | Nonhydrostatic | 4096x1024 (2D) | TEOS10, RK3, ScalarDiffusivity | minutes |
| 6 | Eady LES | Nonhydrostatic | 1024x1024x64 | BackgroundFields, WENO(9) | ~1 hour (1 day) |
| 7 | Headland | Nonhydrostatic | 384x128x64 | ImmersedBoundary, open BCs, TEOS10 | minutes |
| 8 | Internal tides | Hydrostatic | 2000x200 | ImmersedBoundary, WENO, tidal forcing | minutes |
| 9 | Vertical mixing | Hydrostatic | 1x1x50 | CATKE, k-epsilon | seconds |
| 10 | Global baroclinic | Hydrostatic | 1440x680x10 | 3 grid types, WENO-VI, split-explicit | ~1 hour |
| 11 | Coupled ocean-ice | Hydrostatic | 2160x1080x60 | TripolarGrid, ClimaOcean, sea ice | ~hours |

---

## Implementation notes

### Warm-up protocol
Every benchmark script should:
1. Build the model and set initial conditions.
2. Run 10 warm-up time steps (to trigger Julia JIT compilation and CUDA
   kernel compilation).
3. Call `CUDA.@sync` and `GC.gc(true)` to ensure clean state.
4. Record start time with `time_ns()`.
5. Run the timed segment.
6. Call `CUDA.@sync` and record end time.

### Reporting format
Each benchmark should produce a single JSON file with:
```json
{
  "listing": 6,
  "label": "Eady LES",
  "grid_size": [1024, 1024, 64],
  "n_steps_timed": 1440,
  "wall_clock_seconds": 3456.7,
  "time_per_step_ms": 2400.5,
  "sypd": 0.012,
  "peak_gpu_memory_gb": 12.3,
  "gpu": "H100-SXM-80GB",
  "julia_version": "1.11.x",
  "oceananigans_version": "0.9x.y"
}
```

### Dependencies
- Listings 1-8: `Oceananigans`, `CUDA`, `SeawaterPolynomials` (for 5, 7)
- Listing 9: `Oceananigans`, `CUDA`
- Listing 10: `Oceananigans`, `CUDA`, `SeawaterPolynomials`
- Listing 11: `Oceananigans`, `ClimaOcean`, `ClimaSeaIce`, `CUDA`,
  `OrthogonalSphericalShellGrids`, plus data downloads (ECCO, JRA55, ETOPO1)

### Data requirements
- Listings 1-9: No external data needed
- Listing 10: No external data needed (idealized initial conditions)
- Listing 11: Requires ECCO state estimate (~2 GB), JRA55 atmosphere (~50 GB
  for multi-year), ETOPO1 bathymetry (~0.5 GB). Download before benchmarking.

### Order of execution
Run benchmarks in order 1-11. Listings 1-3 serve as validation that the setup
is working correctly. Listings 4-8 benchmark nonhydrostatic and hydrostatic
models at increasing complexity. Listings 9-11 benchmark the hydrostatic model
with parameterizations and coupling at global scale.
