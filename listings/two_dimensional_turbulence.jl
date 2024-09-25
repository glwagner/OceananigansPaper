#=
using Oceananigans
using CairoMakie

grid = RectilinearGrid(CPU(),
                       size=(256, 256), x=(-π, π), y=(-π, π),
                       topology = (Periodic, Periodic, Flat))

function circling_source(x, y, t)
    δ, ω, R = 0.1, 2π/3, 2
    dx = x + R * cos(ω * t)
    dy = y + R * sin(ω * t)
    return exp(-(dx^2 + dy^2) / 2δ^2)
end
                       
model = NonhydrostaticModel(; grid, advection=WENO(order=9),
                            tracers=:c, forcing=(; c=circling_source))

                       
# Set a random initial condition and then run the simulation. 
ϵ(x, y) = 2rand() - 1
set!(model, u=ϵ, v=ϵ)

simulation = Simulation(model; Δt=0.01, stop_time=2.5)

using Printf

function print_progress(sim)
    u, v, w = sim.model.velocities
    max_u = max(maximum(abs, u), maximum(abs, v))
    @info @sprintf("Iteration: %d, time: %.1f, max|u|: %.2e",
                   iteration(sim), time(sim), max_u) 
    return nothing
end

add_callback!(simulation, print_progress, IterationInterval(10))

run!(simulation)
c = deepcopy(model.tracers.c)
simulation.stop_time = 10
run!(simulation)
=#

# Visualize the resulting vorticity field.
u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u))
compute!(ζ)

fig = Figure(size=(800, 400))
axz = Axis(fig[1, 1], aspect=1)
axc = Axis(fig[1, 2], aspect=1)

ζlim = 3/4 * maximum(abs, ζ)
heatmap!(axz, ζ, colormap=:balance, colorrange=(-ζlim, ζlim))
heatmap!(axc, c, colormap=:magma, colorrange=(0, 0.05))

hidedecorations!(axz)
hidedecorations!(axc)

display(fig)

save("first_example.pdf", current_figure())

