# 1° global SYPD benchmark with explicit spinup phase.
# Pattern: spin up at small Δt from cold ECCO IC, then re-time at large Δt.
# Lets us report SYPD at the production Δt (e.g. 30 min) rather than the
# conservative Δt the cold ECCO state needs to avoid NaN.
#
# Env vars:
#   GLOBAL_NX, GLOBAL_NY, GLOBAL_NZ   default 1440, 560, 40
#   GLOBAL_SPINUP_DT_MIN   default 5    (Δt during spinup, in minutes)
#   GLOBAL_SPINUP_DAYS     default 1    (sim-days of spinup)
#   GLOBAL_BENCH_DT_MIN    default 30   (Δt during the timed measurement)

using Oceananigans
using Oceananigans.Units
import NumericalEarth
using Dates
using CFTime
using CUDA
using Printf

include(joinpath(@__DIR__, "sypd_harness.jl"))

const Nx = parse(Int, get(ENV, "GLOBAL_NX", "1440"))
const Ny = parse(Int, get(ENV, "GLOBAL_NY", "560"))
const Nz = parse(Int, get(ENV, "GLOBAL_NZ", "40"))
const spinup_dt_min  = parse(Float64, get(ENV, "GLOBAL_SPINUP_DT_MIN",  "5"))
const spinup_days    = parse(Float64, get(ENV, "GLOBAL_SPINUP_DAYS",    "1"))
const bench_dt_min   = parse(Float64, get(ENV, "GLOBAL_BENCH_DT_MIN",  "30"))

@info @sprintf("Global 1° spinup-then-bench: %d×%d×%d, spinup %g days at Δt=%g min, then bench at Δt=%g min",
               Nx, Ny, Nz, spinup_days, spinup_dt_min, bench_dt_min)

exponential_z_faces(; Nz, depth) = -depth * (1 .- (0:Nz) / Nz).^2
z = exponential_z_faces(; Nz, depth=5000)

arch = GPU()
grid = LatitudeLongitudeGrid(arch; z,
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             longitude = (0, 360),
                             latitude = (-70, 70))

bathymetry = NumericalEarth.regrid_bathymetry(grid)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry))

ocean = NumericalEarth.ocean_simulation(grid)
dates = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T = NumericalEarth.ECCOMetadatum(:temperature; date=dates),
                  S = NumericalEarth.ECCOMetadatum(:salinity; date=dates))

atmosphere = NumericalEarth.JRA55PrescribedAtmosphere(arch)
sea_ice    = NumericalEarth.sea_ice_simulation(grid, ocean)
coupled    = NumericalEarth.OceanSeaIceModel(ocean, sea_ice; atmosphere)

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
if !isfinite(maxT)
    error("[spinup] NaN/Inf in T after spinup — Δt=$(spinup_dt_min) min unstable too. " *
          "Try smaller GLOBAL_SPINUP_DT_MIN.")
end
@info @sprintf("[spinup] post-spinup max|T| = %.3f", maxT)

# === Bench phase: switch to larger Δt and measure ===
simulation.Δt = bench_dt_min * minutes
@info @sprintf("[bench] switching Δt to %g min", bench_dt_min)

measure_sypd!(simulation; name="global_1deg_spinup",
              extra_config=(; Nx, Ny, Nz,
                            spinup_dt_min, spinup_days, bench_dt_min,
                            grid_kind="LatitudeLongitudeGrid"))
