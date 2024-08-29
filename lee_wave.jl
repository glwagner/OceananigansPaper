using Oceananigans
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Units

Nx = 200
Nz = 100
Lx = 200kilometers
Lz = 1000
N² = 1e-4
U₀ = 1

# Set up a grid with a Gaussian mountain bathymetry
underlying_grid = RectilinearGrid(size = (Nx, 1, Nz), halo = (4, 1, 4),
                                  x = (0, Lx), z = (-Lz, 0), y=(0, 1),
                                  topology = (Periodic, Periodic, Bounded))

x₀, h, Δ = Lx/4, Lz/2, Lx/20
bottom_height(x, y) = - Lz + h * exp(-(x - x₀)^2 / 2Δ^2)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

right_side_mask = GaussianMask{:x}(center=Lx, width=Lx/10)
sponge = Relaxation(rate=1/10minutes, mask=right_side_mask)

model = HydrostaticFreeSurfaceModel(; grid, tracers=(:b, :e), buoyancy=BuoyancyTracer(),
                                    closure = CATKEVerticalDiffusivity(),
                                    forcing = (; u=sponge, v=sponge),
                                    momentum_advection=WENO(), tracer_advection=WENO())                                    

bᵢ(x, y, z) = N² * z
set!(model, b=bᵢ, u=U₀)

simulation = Simulation(model; Δt=30, stop_iteration=1000)
run!(simulation)

using CairoMakie

heatmap(model.velocities.u, nan_color=:lightgray, axis=(; aspect=5))
display(current_figure())

