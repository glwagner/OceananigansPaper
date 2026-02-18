# Benchmark Plan for Oceananigans Paper Listings

**Purpose:** Provide wall-clock timing and throughput numbers for a few key examples
in the paper, in response to Reviewer #4 (comment #473) requesting performance
data (wall-clock duration / SYPD) for the code listings.

**Baseline hardware:** Single NVIDIA H100 80 GB SXM GPU.

**General methodology:**

1. Each benchmark script should be derived from the corresponding paper listing
   or repository script, with minimal modifications (shorter runs, disabled I/O).
2. Exclude Julia compilation time by running a short warm-up phase (a few time
   steps) before starting the timing loop.
3. Disable all output writers during the timed portion to isolate compute cost.
4. Report timing as the median of 3 independent runs.
5. Collect GPU memory usage via `CUDA.memory_status()` at peak.

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
- Disable all JLD2 output writers
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

## Summary table

| # | Listing | Model | Grid size | Key features | Expected run time |
|---|---------|-------|-----------|-------------|-------------------|
| 6 | Eady LES | Nonhydrostatic | 1024x1024x64 | BackgroundFields, WENO(9) | ~1 hour (1 day) |
| 7 | Headland | Nonhydrostatic | 384x128x64 | ImmersedBoundary, open BCs, TEOS10 | minutes |
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
- Listings 6-7: `Oceananigans`, `CUDA`, `SeawaterPolynomials` (for 7)
- Listing 10: `Oceananigans`, `CUDA`, `SeawaterPolynomials`
- Listing 11: `Oceananigans`, `ClimaOcean`, `ClimaSeaIce`, `CUDA`,
  `OrthogonalSphericalShellGrids`, plus data downloads (ECCO, JRA55, ETOPO1)

### Data requirements
- Listings 6, 7, 10: No external data needed
- Listing 11: Requires ECCO state estimate (~2 GB), JRA55 atmosphere (~50 GB
  for multi-year), ETOPO1 bathymetry (~0.5 GB). Download before benchmarking.
