using Oceananigans
using Oceananigans.Units

H  = 2kilometers  # domain depth (m)
h₀ = 100          # mountain height (m)
δ  = 20kilometers # mountain width (m)

grid = RectilinearGrid(size = (2000, 200),
                       x = (-1000kilometers, 1000kilometers),
                       z = (-H, 0),
                       halo = (4, 4),
                       topology = (Periodic, Flat, Bounded))

using Random: seed!
seed!(123)
seamounts = 42
W = grid.Lx - 4δ
x₀ = W .* (rand(seamounts) .- 1/2) # mountains' positions ∈ [-Lx/2+2δ, Lx/2-2δ]
h  = h₀ .* (1 .+ rand(seamounts))  # mountains' heights ∈ [h₀, 2h₀]

bottom(x) = -H + sum(h[s] * exp(-(x - x₀[s])^2 / 2δ^2) for s = 1:seamounts)

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))

@inline tidal_forcing(x, z, t, p) = p.U₂ * 2π / p.T₂ * sin(2π / p.T₂ * t)

T₂ = 12.421hours # period of M₂ tide constituent
u_forcing = Forcing(tidal_forcing, parameters=(; U₂=0.1, T₂=T₂))

model = HydrostaticFreeSurfaceModel(; grid, tracers=:b, buoyancy=BuoyancyTracer(),
                                      momentum_advection = WENO(),
                                      tracer_advection = WENO(),
                                      forcing = (; u = u_forcing))

bᵢ(x, z) = 1e-5 * z
set!(model, b=bᵢ)

simulation = Simulation(model; Δt=1minutes, stop_time=16T₂)

writer = JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                          filename = "random_topography_internal_tide.jld2",
                          schedule = TimeInterval(T₂/20),
                          overwrite_existing=true)

simulation.output_writers[:jld2] = writer

run!(simulation)
