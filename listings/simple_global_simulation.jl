using Oceananigans
using Oceananigans.Units
import ClimaOcean
using Dates
using CFTime

Nz = 40
z = ClimaOcean.exponential_z_faces(; Nz, depth=5000)
arch = GPU()
grid = LatitudeLongitudeGrid(arch; z,
                             size = (1440, 560, Nz),
                             halo = (7, 7, 7),
                             longitude = (0, 360),
                             latitude = (-70, 70))

bathymetry = ClimaOcean.regrid_bathymetry(grid) # builds gridded bathymetry based on ETOPO1
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry))

# Build an ocean simulation initialized to the ECCO state estimate on Jan 1, 1993
ocean = ClimaOcean.ocean_simulation(grid)
dates = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T = ClimaOcean.ECCOMetadata(:temperature; dates),
                  S = ClimaOcean.ECCOMetadata(:salinity; dates))

# Build and run an OceanSeaIceModel (with no sea ice component) forced by JRA55 reanalysis
atmosphere = ClimaOcean.JRA55_prescribed_atmosphere(arch)
coupled_model = ClimaOcean.OceanSeaIceModel(ocean; atmosphere)
simulation = Simulation(coupled_model, Δt=5minutes, stop_time=4 * 360days)

using Printf

wall_time = Ref(time_ns())
sim_time = Ref(time(simulation))

function progress(sim)
    u, v, w = sim.model.ocean.model.velocities
    elapsed = 1e-9 * (time_ns() - wall_time[])
    simulated = time(sim) - sim_time[]
    SYPD = simulated / elapsed / 365

    msg = @sprintf("Iter: %d, time: %s, Δt: %s, elapsed: %s, SDPD: %.2f, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(elapsed), SYPD,
                   maximum(abs, u), maximum(abs, v), maximum(abs, w))
    @info msg
    wall_time[] = time_ns()
    sim_time[] = time(sim)
    return nothing
end
                   
add_callback!(simulation, progress, IterationInterval(100))

u, v, w = ocean.model.velocities
ζ = ∂x(v) - ∂y(u)
s = @at (Center, Center, Center) sqrt(u^2 + v^2)
outputs = merge((; ζ, s, w), ocean.model.tracers)
writer = JLD2OutputWriter(ocean.model, outputs,
                          filename = "simple_global_simulation.jld2",
                          indices = (:, :, grid.Nz),
                          array_type = Array{Float32},
                          schedule = TimeInterval(1days),
                          overwrite_existing = true)

checkpointer = Checkpointer(ocean.model,
                            prefix = "simple_global_simulation_checkpoint",
                            schedule = TimeInterval(10days),
                            cleanup = true,
                            overwrite_existing = true)

simulation.output_writers[:surface] = writer
simulation.output_writers[:chk] = checkpointer

run!(simulation)

