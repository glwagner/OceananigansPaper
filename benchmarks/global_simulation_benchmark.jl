# Global ocean-sea-ice SYPD benchmark.
# Mirrors listings/simple_global_simulation.jl but uses the SYPD harness for
# a short timed measurement window. The listing config is 1° lat-lon
# (1440×560×40); the paper plan calls for 1/6° tripolar (2160×1080×60),
# which is multi-GPU territory — that's a separate script.
#
# Requires ECCO_USERNAME + ECCO_WEBDAV_PASSWORD env vars (NASA Earthdata).
# Pre-stage data with: julia --project scripts/prestage_data.jl

using Oceananigans
using Oceananigans.Units
import NumericalEarth
using Dates
using CFTime
using CUDA
using Printf

include(joinpath(@__DIR__, "sypd_harness.jl"))

const Nz = parse(Int, get(ENV, "GLOBAL_NZ", "40"))
const Nx = parse(Int, get(ENV, "GLOBAL_NX", "1440"))
const Ny = parse(Int, get(ENV, "GLOBAL_NY", "560"))
const dt_minutes = parse(Float64, get(ENV, "GLOBAL_DT_MIN", "5"))

@info @sprintf("Global benchmark config: %d×%d×%d at Δt=%g min", Nx, Ny, Nz, dt_minutes)

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
sea_ice = NumericalEarth.sea_ice_simulation(grid, ocean)
coupled_model = NumericalEarth.OceanSeaIceModel(ocean, sea_ice; atmosphere)

simulation = Simulation(coupled_model, Δt = dt_minutes * minutes)

measure_sypd!(simulation; name="simple_global_simulation",
              extra_config=(; Nx, Ny, Nz, dt_minutes,
                            grid_kind="LatitudeLongitudeGrid",
                            atmosphere="JRA55", ocean_state="ECCO 1993-01-01"))
