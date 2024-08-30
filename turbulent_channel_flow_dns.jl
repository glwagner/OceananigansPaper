using Oceananigans

grid = RectilinearGrid(GPU(),
                       size=(512, 256, 256),
                       x = (0, 4π)
                       y = (-π, π)
                       z = (-1, 1)
                       topology = (Periodic, Periodic, Bounded))

no_slip = ValueBoundaryCondition(0)
u_bcs = FieldBoundaryConditions(top=no_slip, bottom=no_slip)
v_bcs = FieldBoundaryConditions(top=no_slip, bottom=no_slip)
boundary_conditions = (u=u_bcs, v=v_bcs)

Re = 180
forcing = (; u=1)
closure = ScalarDiffusivity(ν=1/Re)
model = NonhydrostaticModel(; grid, closure, boundary_conditions, forcing)

u₀(x, y, z) = 30 * (1 - z^2)
v₀(x, y, z) = rand() * 0.001
w₀(x, y, z) = rand() * 0.001

set!(model, u = u₀, v = v₀, w = w₀)

simulation = Simulation(model, Δt=1e-4, stop_time=5)

run!(simulation)

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




