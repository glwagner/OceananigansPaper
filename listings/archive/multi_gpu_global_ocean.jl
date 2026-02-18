using ClimaOcean
using Oceananigans
using Oceananigans.Units
using CFTime
using Dates
using Printf

arch = Distributed(GPU(), partition=Partition(8), synchronized_communication=true)

Nx = 4320
Ny = 1800
Nz = 40

depth = 6000meters
z = ExponentialDiscretization(Nz, -depth, 0)

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z,
                             latitude  = (-75, 75),
                             longitude = (0, 360))

bottom_height = regrid_bathymetry(grid; 
                                  minimum_depth = 10meters,
                                  interpolation_passes = 5,
                                  major_basins = 3)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true)

momentum_advection = WENOVectorInvariant()
tracer_advection = WENO(order=7)

free_surface = SplitExplicitFreeSurface(grid; substeps = 70)
ocean = ocean_simulation(grid; free_surface, momentum_advection, tracer_advection)

date = DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T=ECCOMetadata(:temperature; dates=date),
                  S=ECCOMetadata(:salinity; dates=date))

radiation = Radiation(arch)
atmosphere    = JRA55PrescribedAtmosphere(arch; backend=JRA55NetCDFBackend(41))
coupled_model = OceanSeaIceModel(ocean; atmosphere, radiation)

simulation = Simulation(coupled_model; Δt=30, stop_time=30days)

wall_time = Ref(time_ns())

function progress(sim) 
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities
    T = ocean.model.tracers.T

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))

    umax = (maximum(abs, interior(u)),
            maximum(abs, interior(v)),
            maximum(abs, interior(w)))

    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Iter: %d, time: %s, Δt: %s", iteration(sim), prettytime(sim), prettytime(sim.Δt))
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): (%.2f, %.2f) ᵒC, wall time: %s",
                    umax..., Tmax, Tmin, prettytime(step_time))

    ClimaOcean.@root @info msg 

    wall_time[] = time_ns()
end

simulation.callbacks[:progress] = Callback(progress, TimeInterval(5days))

outputs = merge(ocean.model.tracers, ocean.model.velocities)
ocean.output_writers[:surface] = JLD2Writer(ocean.model, outputs;
                                            schedule = TimeInterval(1days),
				                            filename = "near_global_surface_fields_$(arch.local_rank)",
                                            indices = (:, :, grid.Nz),
                                            with_halos = true,
                                            overwrite_existing = true,
                                            array_type = Array{Float32})

run!(simulation)
simulation.stop_time = 120days
simulation.Δt = 3.5minutes
run!(simulation)
