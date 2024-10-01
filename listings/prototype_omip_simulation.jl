using Oceananigans
using Oceananigans.Units

using ClimaOcean
using ClimaOcean.ECCO: ECCO_restoring_forcing, ECCO4Monthly, ECCO2Daily, ECCOMetadata

using OrthogonalSphericalShellGrids
using Printf

using CFTime
using Dates

resolution = 1/6 # degree
Nx = round(Int, 360 / resolution)
Ny = round(Int, 170 / resolution)
Nz = 10
Δt = 2minutes
stop_time = 10days
z_faces = exponential_z_faces(; Nz, depth=6000)
arch = GPU() #Distributed(GPU(), partition = Partition(2))

grid = TripolarGrid(arch; 
                    size = (Nx, Ny, Nz), 
                    halo = (7, 7, 7), 
                    z = z_faces, 
                    north_poles_latitude = 55,
                    first_pole_longitude = 75)

bottom_height = retrieve_bathymetry(grid;
                                    minimum_depth = 10,
                                    interpolation_passes = 2,
                                    major_basins = 1)
 
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height); active_cells_map=true) 
free_surface = SplitExplicitFreeSurface(grid; cfl=0.7, fixed_Δt=Δt)

#####
##### Add restoring to ECCO fields for temperature and salinity in the artic and antarctic
#####

@inline function restoring_mask(λ, φ, z, t=0)
    ϵN = (φ - 75) / 5
    ϵN = clamp(ϵN, zero(ϵN), one(ϵN))
    ϵS = - (φ + 75) / 5
    ϵS = clamp(ϵS, zero(ϵS), one(ϵS))
    return ϵN + ϵS
end

restoring_mask_field = CenterField(grid)
set!(restoring_mask_field, restoring_mask)

@inline sponge_layer(λ, φ, z, t, c, ω) = - restoring_mask(λ, φ, z, t) * ω * c
Fu = Forcing(sponge_layer, field_dependencies=:u, parameters=1/1days)
Fv = Forcing(sponge_layer, field_dependencies=:v, parameters=1/1days)

dates = DateTimeProlepticGregorian(1993, 1, 1) : Month(1) : DateTimeProlepticGregorian(2003, 12, 1)

temperature = ECCOMetadata(:temperature, dates, ECCO4Monthly())
salinity    = ECCOMetadata(:salinity,    dates, ECCO4Monthly())

FT = ECCO_restoring_forcing(temperature; mask=restoring_mask_field, grid, architecture=arch, timescale=1days)
FS = ECCO_restoring_forcing(salinity;    mask=restoring_mask_field, grid, architecture=arch, timescale=1days)
forcing = (; T=FT, S=FS, u=Fu, v=Fv)

# New advection scheme
tracer_advection = WENO(order=5)
momentum_advection = WENOVectorInvariant(vorticity_order=5)

ocean = ocean_simulation(grid; forcing, tracer_advection, free_surface, momentum_advection) 

backend = JRA55NetCDFBackend(24) 
atmosphere = JRA55_prescribed_atmosphere(arch; backend)
radiation = Radiation(arch, ocean_albedo=LatitudeDependentAlbedo())
sea_ice = ClimaOcean.OceanSeaIceModels.MinimumTemperatureSeaIce()
coupled_model = OceanSeaIceModel(ocean, sea_ice; atmosphere, radiation)

@info "Set up coupled model:"
@info coupled_model

coupled_simulation = Simulation(coupled_model; Δt, stop_time)

prefix = "prototype_omip_simulation"

#=
fluxes = (u = ocean.model.velocities.u.boundary_conditions.top.condition,
          v = ocean.model.velocities.v.boundary_conditions.top.condition,
          T = ocean.model.tracers.T.boundary_conditions.top.condition,
          S = ocean.model.tracers.S.boundary_conditions.top.condition)

ocean.output_writers[:fluxes] = JLD2OutputWriter(ocean.model, fluxes,
                                                 schedule = TimeInterval(1days),
                                                 with_halos = true,
                                                 overwrite_existing = true,
                                                 array_type = Array{Float32},
                                                 filename = prefix * "surface_fluxes")
=#

outputs = merge(ocean.model.tracers, ocean.model.velocities)

ocean.output_writers[:surface] = JLD2OutputWriter(ocean.model, outputs, 
                                                  schedule = TimeInterval(1days),
                                                  with_halos = true,
                                                  overwrite_existing = true,
                                                  array_type = Array{Float32},
                                                  filename = prefix * "_surface",
                                                  indices = (:, :, grid.Nz))

#=
ocean.output_writers[:snapshots] = JLD2OutputWriter(ocean.model, merge(ocean.model.tracers, ocean.model.velocities),
                                                    schedule = TimeInterval(10days),
                                                    with_halos = true,
                                                    overwrite_existing = true,
                                                    array_type = Array{Float32},
                                                    filename = prefix * "snapshots")

ocean.output_writers[:checkpoint] = Checkpointer(ocean.model, 
                                                 schedule = TimeInterval(60days),
                                                 overwrite_existing = true,
                                                 prefix = prefix * "checkpoint")
=#

@info "Built an ocean simulation with the model:"
@info ocean.model
@info "\ninside the simulation:"
@info ocean

#####
##### The atmosphere
#####

# Spin up
set!(ocean.model, 
     T = ECCOMetadata(:temperature, dates[1], ECCO2Daily()),
     S = ECCOMetadata(:salinity,    dates[1], ECCO2Daily()))

wall_time = Ref(time_ns())

function progress(sim) 
    ocean = sim.model.ocean
    u, v, w = ocean.model.velocities  
    T, S = ocean.model.tracers

    Tmax = maximum(interior(T))
    Tmin = minimum(interior(T))
    umax = maximum(abs, interior(u)), maximum(abs, interior(v)), maximum(abs, interior(w))
    step_time = 1e-9 * (time_ns() - wall_time[])

    msg = @sprintf("Time: %s, step: %d", prettytime(sim), iteration(sim))

    msg *= @sprintf(", max(u): (%.2e, %.2e, %.2e) m s⁻¹, extrema(T): %.2fᵒC, %.2fᵒC, wall time: %s",
                   umax..., Tmax, Tmin, prettytime(step_time))

    @info msg
    wall_time[] = time_ns()
end

add_callback!(coupled_simulation, progress, IterationInterval(10))

@info "Running the coupled simulation:"
@info coupled_simulation
    
run!(coupled_simulation)

# Run for real
#coupled_simulation.Δt = 5minutes
#
## Let's reset the maximum number of iterations
#coupled_simulation.stop_time = 7200days
#coupled_simulation.stop_iteration = Inf
#
#run!(coupled_simulation)

