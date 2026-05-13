# 1/6° tripolar global ocean–sea-ice SYPD benchmark.
# Paper plan config:
#   Grid: Tripolar, 1/6° (2160 × 1080 × 60), exponential z to 6 km
#   Momentum advection: WENO vector invariant (order=9)
#   Tracer advection: WENO(order=7)
#   Vertical mixing: CATKE
#   EOS: TEOS10
#   Sea ice: ClimaSeaIce (sea_ice_simulation default — 0-layer slab + viscoplastic dynamics)
#   Atmosphere: JRA55 prescribed
#   Init: ECCO Jan 1 1993
#   Δt: 6 min (production rate; spinup at 20 s not benchmarked here)
#
# Requires ECCO_USERNAME + ECCO_WEBDAV_PASSWORD. Pre-stage with
# `julia --project scripts/prestage_data.jl` (the 1° stage covers this since
# JRA55 is shared and ETOPO is regridded per-grid).

using Oceananigans
using Oceananigans.Units
using Oceananigans.OrthogonalSphericalShellGrids: TripolarGrid
using Oceananigans.Advection: WENOVectorInvariant
import NumericalEarth
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Dates
using CFTime
using CUDA
using Printf

include(joinpath(@__DIR__, "sypd_harness.jl"))

const Nx = parse(Int, get(ENV, "GLOBAL16_NX", "2160"))
const Ny = parse(Int, get(ENV, "GLOBAL16_NY", "1080"))
const Nz = parse(Int, get(ENV, "GLOBAL16_NZ", "60"))
const dt_minutes = parse(Float64, get(ENV, "GLOBAL16_DT_MIN", "6"))

@info @sprintf("Global 1/6° benchmark config: tripolar %d×%d×%d at Δt=%g min", Nx, Ny, Nz, dt_minutes)

exponential_z_faces(; Nz, depth) = -depth * (1 .- (0:Nz) / Nz).^2
z = exponential_z_faces(; Nz, depth=6000)

arch = GPU()
grid = TripolarGrid(arch; size=(Nx, Ny, Nz), z, halo=(7, 7, 7))

bathymetry = NumericalEarth.regrid_bathymetry(grid)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry))

ocean = NumericalEarth.ocean_simulation(grid;
    momentum_advection = WENOVectorInvariant(order=9),
    tracer_advection   = WENO(order=7),
    equation_of_state  = TEOS10EquationOfState(),
)

dates = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T = NumericalEarth.ECCOMetadatum(:temperature; date=dates),
                  S = NumericalEarth.ECCOMetadatum(:salinity; date=dates))

atmosphere = NumericalEarth.JRA55PrescribedAtmosphere(arch)
sea_ice    = NumericalEarth.sea_ice_simulation(grid, ocean)
coupled    = NumericalEarth.OceanSeaIceModel(ocean, sea_ice; atmosphere)

simulation = Simulation(coupled, Δt = dt_minutes * minutes)

measure_sypd!(simulation; name="global_16th_tripolar",
              extra_config=(; Nx, Ny, Nz, dt_minutes,
                            grid_kind="TripolarGrid", depth_m=6000,
                            momentum_advection="WENOVectorInvariant(9)",
                            tracer_advection="WENO(7)",
                            atmosphere="JRA55", ocean_state="ECCO 1993-01-01"))
