# 1/6° tripolar global SYPD benchmark with explicit spinup phase.
# Spin up at Δt=spinup_dt_min from cold ECCO IC, then re-time at Δt=bench_dt_min.
# Lets us push past the paper-plan Δt=6 min and find the largest stable Δt
# for the WENOVectorInvariant(9) momentum + WENO(7) tracer config.
#
# Env vars:
#   GLOBAL16_SPINUP_DT_MIN   default 6     (Δt during spinup)
#   GLOBAL16_SPINUP_DAYS     default 1     (sim-days of spinup)
#   GLOBAL16_BENCH_DT_MIN    default 10    (Δt during the timed measurement)

using Oceananigans
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: TripolarGrid
using Oceananigans.Advection: WENOVectorInvariant
using NumericalEarth
using NumericalEarth: Metadatum, ECCO4Monthly, JRA55PrescribedAtmosphere,
                      JRA55NetCDFBackend, Radiation, regrid_bathymetry,
                      ocean_simulation, sea_ice_simulation, OceanSeaIceModel
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Dates
using CFTime
using CUDA
using Printf

include(joinpath(@__DIR__, "sypd_harness.jl"))

const Nx = parse(Int, get(ENV, "GLOBAL16_NX", "2160"))
const Ny = parse(Int, get(ENV, "GLOBAL16_NY", "1080"))
const Nz = parse(Int, get(ENV, "GLOBAL16_NZ", "60"))
const spinup_dt_min  = parse(Float64, get(ENV, "GLOBAL16_SPINUP_DT_MIN", "6"))
const spinup_days    = parse(Float64, get(ENV, "GLOBAL16_SPINUP_DAYS",   "1"))
const bench_dt_min   = parse(Float64, get(ENV, "GLOBAL16_BENCH_DT_MIN", "10"))

@info @sprintf("Global 1/6° spinup-then-bench: %d×%d×%d tripolar, spinup %g days at Δt=%g min, bench at Δt=%g min",
               Nx, Ny, Nz, spinup_days, spinup_dt_min, bench_dt_min)

arch  = GPU()
depth = 6000meters
z = ExponentialDiscretization(Nz, -depth, 0; mutable=true)   # z* coord

raw_grid = TripolarGrid(arch; size=(Nx, Ny, Nz), z, halo=(7, 7, 7))
bottom_height = regrid_bathymetry(raw_grid;
    minimum_depth=10meters, interpolation_passes=5, major_basins=3)
grid = ImmersedBoundaryGrid(raw_grid, GridFittedBottom(bottom_height); active_cells_map=true)

# substeps=100 keeps the substep COUNT fixed so the substep size shrinks with
# Δt; this means an in-place `simulation.Δt = bench_dt` is safe — no model
# rebuild needed and the gravity-wave CFL stays satisfied at both spinup_dt
# (= 5 min, 3 s/substep) and bench_dt (= 15-20 min, 9-12 s/substep).
const fs_substeps = parse(Int, get(ENV, "GLOBAL16_SUBSTEPS", "100"))
free_surface = SplitExplicitFreeSurface(grid; substeps=fs_substeps)
@info @sprintf("free_surface: SplitExplicitFreeSurface(substeps=%d)", fs_substeps)

ocean = ocean_simulation(grid;
    momentum_advection = WENOVectorInvariant(order=9),
    tracer_advection   = WENO(order=7),
    equation_of_state  = TEOS10EquationOfState(),
    free_surface       = free_surface,
)

set!(ocean.model, T = Metadatum(:temperature, dataset=ECCO4Monthly()),
                  S = Metadatum(:salinity,    dataset=ECCO4Monthly()))

atmosphere = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))
radiation  = Radiation(arch)
sea_ice    = sea_ice_simulation(grid, ocean; advection=WENO(order=5))
set!(sea_ice.model, h = Metadatum(:sea_ice_thickness,     dataset=ECCO4Monthly()),
                    ℵ = Metadatum(:sea_ice_concentration, dataset=ECCO4Monthly()))

coupled    = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)
simulation = Simulation(coupled, Δt = spinup_dt_min * minutes)

# === Spinup phase ===
spinup_iters = round(Int, spinup_days * 86400 / (spinup_dt_min * 60))
@info @sprintf("[spinup] running %d iters at Δt=%g min (= %g sim-days)",
               spinup_iters, spinup_dt_min, spinup_days)
simulation.stop_iteration = spinup_iters
simulation.stop_time = Inf
empty!(simulation.output_writers)
empty!(simulation.diagnostics)
spinup_wall_start = time()
run!(simulation)
spinup_wall = time() - spinup_wall_start
@info @sprintf("[spinup] done: %d iters in %.1f s wall (%.3f sec/iter)",
               simulation.model.clock.iteration, spinup_wall,
               spinup_wall / max(1, simulation.model.clock.iteration))

# Sanity check before switching Δt
T = simulation.model.ocean.model.tracers.T
maxT = maximum(abs, T)
isfinite(maxT) || error("[spinup] NaN/Inf in T after spinup at Δt=$(spinup_dt_min) min — pick smaller GLOBAL16_SPINUP_DT_MIN.")
@info @sprintf("[spinup] post-spinup max|T| = %.3f", maxT)

# === Bench phase: switch to candidate Δt and measure ===
simulation.Δt = bench_dt_min * minutes
@info @sprintf("[bench] switching Δt to %g min", bench_dt_min)

measure_sypd!(simulation; name="global_16th_spinup",
              extra_config=(; Nx, Ny, Nz,
                            spinup_dt_min, spinup_days, bench_dt_min,
                            grid_kind="TripolarGrid",
                            momentum_advection="WENOVectorInvariant(9)",
                            tracer_advection="WENO(7)"))
