using Oceananigans
using Oceananigans.Units

grid = LatitudeLongitudeGrid(GPU(), size=(32 * 20, 32 * 22, 2),
                             longitude = (0, 22), latitude = (30, 50), z = (-2000, 0))

@inline τx(λ, φ, t, τ₀) = τ₀ * (1 - cos(2π * (φ - 30) / 20))
u_top_bc = FluxBoundaryCondition(τx, parameters=1e-4)
u_bcs = FieldBoundaryConditions(top=u_top_bc)

coriolis = HydrostaticSphericalCoriolis()

model = HydrostaticFreeSurfaceModel(; grid, coriolis, tracers=:b, buoyancy=BuoyancyTracer(),
                                    boundary_conditions = (; u_bcs))

