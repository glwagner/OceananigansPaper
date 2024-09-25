using Oceananigans
using Oceananigans.Units

H  = 2kilometers  # domain depth (m)
h₁ = 500          # first mountain height (m)
h₂ = 1000         # second mountain height (m)
δ  = 20kilometers # mountain width (m)

underlying_grid = RectilinearGrid(size = (500, 100),
                                  x = (-1000kilometers, 1000kilometers),
                                  z = (-H, 0),
                                  halo = (4, 4),
                                  topology = (Periodic, Flat, Bounded))

N² = 1e-5
f = 1e-4
@show λ = π * sqrt(N²) * H / f # first mode wavelength
x₁ = + 2λ# / 2
x₂ = - 2λ# / 2

bottom(x) = -H + h₁ * exp(-(x - x₁)^2 / 2δ^2) +
                 h₂ * exp(-(x - x₂)^2 / 2δ^2)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))

@inline tidal_forcing(x, z, t, p) = p.U₂ * 2π/p.T₂ * sin(2π / p.T₂ * t)
u_forcing = Forcing(tidal_forcing, parameters=(; U₂=0.1, T₂=12.421hours))

model = HydrostaticFreeSurfaceModel(; grid, tracers=:b, buoyancy=BuoyancyTracer(),
                                    coriolis = FPlane(; f),
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO(),
                                    forcing = (; u = u_forcing))

bᵢ(x, z) = 1e-5 * z
set!(model, b=bᵢ)
T₂ = 12.421hours
simulation = Simulation(model; Δt=2minutes, stop_time=8T₂)
run!(simulation)

using GLMakie

fig = Figure(size=(1200, 400))
ax = Axis(fig[1, 1], xlabel="x (km)", ylabel="z (m)", aspect=6)
x, y, z = nodes(model.velocities.w)
heatmap!(ax, x ./ 1e3, z, model.velocities.w, colormap=:balance, nan_color=:gray)

