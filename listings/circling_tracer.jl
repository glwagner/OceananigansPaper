using Oceananigans

# A two-dimensional square grid on the CPU.
grid = RectilinearGrid(CPU(),
                       size = (256, 256),
                       x = (-π, π),
                       y = (-π, π),
                       topology = (Periodic, Periodic, Flat))

function circling_source(x, y, t)
    δ, ω, R = 0.05, 2π/3, 2
    dx = x + R * cos(ω * t)
    dy = y + R * sin(ω * t)
    return exp(-(dx^2 + dy^2) / 2δ^2)
end
                       
# Build a model with the "Weighted Essentially Non-Oscillatory" (WENO)
# advection scheme, which introduces implicit numerical diffusion.
model = NonhydrostaticModel(; grid, advection=WENO(),
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

using CairoMakie

axis = (; aspect=1, title="Passive tracer at t=2.5", xlabel="x", ylabel="y")
heatmap(model.tracers.c, colormap=:thermal; axis)

save("second_example.pdf", current_figure())

