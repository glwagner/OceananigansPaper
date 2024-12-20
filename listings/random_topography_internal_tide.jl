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
positions = W .* (rand(seamounts) .- 1/2) # ∈ [-Lx, Lx]
heights = h₀ .* (1 .+ rand(seamounts)) # ∈ [h₀, 2h₀]

function bottom(x)
    z = - H
    for (h, ξ) in zip(heights, positions)
        z += h * exp(-(x - ξ)^2 / 2δ^2)
    end
    return z
end

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))

@inline tidal_forcing(x, z, t, p) = p.U₂ * 2π/p.T₂ * sin(2π / p.T₂ * t)
u_forcing = Forcing(tidal_forcing, parameters=(; U₂=0.1, T₂=12.421hours))

model = HydrostaticFreeSurfaceModel(; grid, tracers=:b, buoyancy=BuoyancyTracer(),
                                      momentum_advection = WENO(),
                                      tracer_advection = WENO(),
                                      forcing = (; u = u_forcing))

bᵢ(x, z) = 1e-5 * z
set!(model, b=bᵢ)
T₂ = 12.421hours
simulation = Simulation(model; Δt=1minutes, stop_time=16T₂)

writer = JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                          filename = "random_topography_internal_tide.jld2",
                          schedule = TimeInterval(T₂/20),
                          overwrite_existing=true)

simulation.output_writers[:jld2] = writer

run!(simulation)

