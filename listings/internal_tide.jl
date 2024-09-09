using Oceananigans
using Oceananigans.Units

H = 2kilometers  # domain depth (m)
h = 500          # mountain height (m)
δ = 20kilometers # mountain width (m)
N² = 1e-4        # initial buoyancy frequency (s⁻²)

underlying_grid = RectilinearGrid(size = (500, 200),
                                  x = (-1000kilometers, 1000kilometers),
                                  z = (-H, 0),
                                  halo = (4, 4),
                                  topology = (Periodic, Flat, Bounded))

bottom(x) = - H + h * exp(-x^2 / 2δ^2)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))

@inline tidal_forcing(x, z, t, p) = p.U₂ * 2π/p.T₂ * sin(2π / p.T₂ * t)
u_forcing = Forcing(tidal_forcing, parameters=(; U₂=0.2, T₂=12.421hours))

model = HydrostaticFreeSurfaceModel(; grid, tracers=:b, buoyancy=BuoyancyTracer(),
                                      momentum_advection = WENO(),
                                      tracer_advection = WENO(),
                                      forcing = (; u = u_forcing))

bᵢ(x, z) = N² * z
set!(model, b=bᵢ)
simulation = Simulation(model; Δt=1minutes, stop_time=2days)
run!(simulation)

using GLMakie

fig = Figure(size=(1200, 400))
ax = Axis(fig[1, 1], xlabel="x (km)", ylabel="z (m)")
x, y, z = nodes(model.velocities.w)
heatmap!(ax, x ./ 1e3, z, model.velocities.w, colormap=:balance, nan_color=:gray)

