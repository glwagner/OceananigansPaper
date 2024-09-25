using Oceananigans
using Oceananigans.Units
using GLMakie

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

coriolis = HydrostaticSphericalCoriolis()

model = HydrostaticFreeSurfaceModel(; grid, coriolis,
                                    momentum_advection = WENOVectorInvariant(),
                                    tracer_advection = WENO(order=9),
                                    tracers = :b,
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; u=u_bcs))

N² = 1e-5
bᵢ(x, y, z) = N² * z
set!(model, b=bᵢ)

simulation = Simulation(model, Δt=10minutes, stop_time=30days)

progress(sim) = @info string("Iter: ", iteration(sim), ", time: ", prettytime(sim))
add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)

Nz = size(grid, 3)
heatmap(view(model.velocities.v, :, :, Nz)) #grid.immersed_boundary.bottom_height)

