using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Printf
# using GLMakie

arch = GPU()
Nx = 8 * 22
Ny = 8 * 20
Nz = 16

grid = LatitudeLongitudeGrid(arch,
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             longitude = (0, 22),
                             latitude = (30, 50),
                             z = (-2000, 0))

δλ = 10
δφ = 20
paraboloid(λ, φ) = (λ - 11)^2 / δλ^2 + (φ - 30)^2 / δφ^2

H∞ = 3000
bottom_height(λ, φ) = - H∞ * (1.1 - paraboloid(λ, φ))
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

@inline τx(λ, φ, t, τ₀) = -τ₀ * (1 - cos(2π * (φ - 30) / 20))
u_top_bc = FluxBoundaryCondition(τx, parameters=1e-4)
u_bcs = FieldBoundaryConditions(top=u_top_bc)

model = HydrostaticFreeSurfaceModel(; grid,
                                    coriolis = HydrostaticSphericalCoriolis(),
                                    momentum_advection = WENOVectorInvariant(),
                                    tracer_advection = WENO(order=9),
                                    closure = CATKEVerticalDiffusivity(),
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; u=u_bcs))

N² = 1e-5
bᵢ(x, y, z) = N² * z
set!(model, b=bᵢ)

simulation = Simulation(model, Δt=1minutes, stop_time=30days)

function progress(sim)
    u, v, w = sim.model.velocities
    @info @sprintf("Iter: %d, time: %s, max|w|: %.2e",
                   iteration(sim), prettytime(sim), maximum(w))

    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

e = model.tracers.e
b = model.tracers.b
u, v, w = model.velocities
κc = model.diffusivity_fields.κc
outputs = (; u, v, w, b, e, κc)

Nx, Ny, Nz = size(grid)
i = floor(Int, Nx/2)
j = floor(Int, Ny/2)
indiceses = [(:, j, :), (i, :, :),
             (:, :, Nz), (:, :, 1), (:, :, Nz-8)]
names = [:xz, :yz, :xy1, :xy2, :xy3]
prefix = "double_gyre_Nx$(Nx)_Nz$(Nz)"

for (name, indices) in zip(names, indiceses)
    output_writer = JLD2OutputWriter(model, outputs; indices,
                                     filename = string(prefix, "_", name),
                                     schedule = TimeInterval(3hours),
                                     with_halos = true,
                                     overwrite_existing = true)
                                          
    simulation.output_writers[name] = output_writer
end

run!(simulation)

#Nz = size(grid, 3)
# heatmap(view(model.velocities.v, :, :, Nz)) #grid.immersed_boundary.bottom_height)
# 
