using Oceananigans
using Oceananigans.Fields: ConstantField
using Oceananigans.TurbulenceClosures:  VerticallyImplicitTimeDiscretization

using Printf

Ny = 256
Nx = 2Ny
Nz = Ny

ϵ = 1 # stretching parameter (↑↑ ⟹  more stretching)
ζ(k) = 2/Nz * (k - 1) - 1 # ζ ∈ [-1, 1]
z(k) = tanh(ϵ * ζ(k)) / tanh(ϵ)

grid = RectilinearGrid(GPU(); topology = (Periodic, Periodic, Bounded),
                       size=(Nx, Ny, Nz), x = (0, 4π), y = (-π, π), z)

@show grid

vitd = VerticallyImplicitTimeDiscretization()
Re = 180
closure = ScalarDiffusivity(vitd; ν=1/Re)

Px = ConstantField(1) # mean x-direction pressure gradient
forcing = (; u=Px)

no_slip = ValueBoundaryCondition(0)
uv_bcs = FieldBoundaryConditions(top=no_slip, bottom=no_slip)
boundary_conditions = (u=uv_bcs, v=uv_bcs)

model = NonhydrostaticModel(; grid, closure, boundary_conditions, forcing, timestepper=:RungeKutta3)

u₀(x, y, z) = 10 * (1 - z^2) + randn()
v₀(x, y, z) = randn()
w₀(x, y, z) = randn()
set!(model, u=u₀, v=v₀, w=w₀)

Δx = minimum_xspacing(grid)
max_Δt = 0.2 * Δx^2 * Re
simulation = Simulation(model, Δt=1e-4, stop_time=10)
conjure_time_step_wizard!(simulation, cfl=0.5, IterationInterval(10); max_Δt)

progress(sim) = @info @sprintf("Iter: %d, time: %.4f, Δt: %.2e, max|u|: (%.2e, %.2e, %.2e)",
                               iteration(sim), time(sim), sim.Δt,
                               maximum(abs, interior(sim.model.velocities.u)),
                               maximum(abs, interior(sim.model.velocities.v)),
                               maximum(abs, interior(sim.model.velocities.w)))

add_callback!(simulation, progress, IterationInterval(100))

simulation.output_writers[:jld2] = JLD2OutputWriter(model, model.velocities;
                                                    filename = "turbulent_channel_dns.jld2",
                                                    schedule = TimeInterval(10),
                                                    with_halos = true,
                                                    overwrite_existing = true)

run!(simulation)

#=
# Probably we need to average in time.
# If results are ok with just space average we can keep them.
u, v, w = model.velocities
U = Field(Average(u, dims=(1, 2)))

u′  = u - U
u′² = Field(u′ * u′)
v′² = Field(v * v)
w′² = Field(w * w)

compute!(u′²)
compute!(v′²)
compute!(w′²)

using CairoMakie

zC⁺ = znodes(u′²)[1:128] .* Re 
zF⁺ = znodes(w′²)[1:128] .* Re 

fig = Figure()
ax  = Axis(fig[1, 1], title = "Mean flow", xscale = log10)
lines!(axm zC⁺, view(U, 1, 1, 1:128))
ax  = Axis(fig[1, 1], title = "Turbulent statistics", xscale = log10) 
lines!(ax, zC⁺, view(u′², 1, 1, 1:128), label = "u′²")
lines!(ax, zC⁺, view(v′², 1, 1, 1:128), label = "v′²")
lines!(ax, zF⁺, view(w′², 1, 1, 1:128), label = "w′²")
=#
