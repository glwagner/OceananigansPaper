using Oceananigans
using ClimaOcean
using OrthogonalSphericalShellGrids
using CFTime
using Dates

Nx = 2160
Ny = 1080
Nz = 60
z = exponential_z_faces(Nz=60, depth=6000)

grid = TripolarGrid(GPU; size = (Nx, Ny, Nz), halo = (7, 7, 7), z)

bottom_height = retrieve_bathymetry(grid, minimum_depth = 5, major_basins = 1,
                                    interpolation_passes = 10)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

start_date = DateTimeProlepticGregorian(1993, 1, 1)
end_date = DateTimeProlepticGregorian(2003, 12, 1)
dates = range(start_date, stop=end_date, step=Month(1))

FT = ECCORestoring(arch, :temperature; dates, rate=1/30days)
FS = ECCORestoring(arch, :salinity; dates, rate=1/30days)
forcing = (; T=FT, S=FS)

ocean = ocean_simulation(grid; free_surface, forcing) 

backend    = JRA55NetCDFBackend(24) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
sea_ice = MinimumTemperatureSeaIce()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere)

# Spin up for a week
simulation = Simulation(coupled_model; Δt=30, stop_time=7days)
run!(simulation)

# Longer run with a long time step
simulation.Δt = 5minutes
simulation.stop_time = 360days
run!(simulation)

