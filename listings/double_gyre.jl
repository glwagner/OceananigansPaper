using Oceananigans
using Oceananigans.Units

Nx = 32 * 22
Ny = 32 * 20
Nz = 2

grid = LatitudeLongitudeGrid(GPU(), size=(Nx, Ny, 2),
                             longitude = (0, 22),
                             latitude = (30, 50),
                             z = (-2000, 0))

δλ = 5
δφ = 20
H∞ = 4000
bottom_height(λ, φ) = -H + (λ - 11)^2 / δλ^2 + φ^2 / δφ^2

@inline τx(λ, φ, t, τ₀) = τ₀ * (1 - cos(2π * (φ - 30) / 20))
u_top_bc = FluxBoundaryCondition(τx, parameters=1e-4)
u_bcs = FieldBoundaryConditions(top=u_top_bc)

coriolis = HydrostaticSphericalCoriolis()

model = HydrostaticFreeSurfaceModel(; grid, coriolis,
                                    tracers = :b,
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (; u_bcs))

