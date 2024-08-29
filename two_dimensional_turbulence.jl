using Oceananigans

# A two-dimensional square grid on the CPU.
grid = RectilinearGrid(CPU(),
                       size = (256, 256),
                       x = (-π, π),
                       y = (-π, π),
                       topology = (Periodic, Periodic, Flat))
                       
# Build a model with the "Weighted Essentially Non-Oscillatory" (WENO)
# advection scheme, which introduces implicit numerical diffusion.
model = NonhydrostaticModel(; grid, advection=WENO())

# Set a random initial condition and then run the simulation. 
ϵ(x, y) = 2rand() - 1
set!(model, u=ϵ, v=ϵ)

simulation = Simulation(model; Δt=0.01, stop_time=10)

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

# Visualize the resulting vorticity field.
u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u))
compute!(ζ)

using CairoMakie

ζlim = 3/4 * maximum(abs, ζ)
axis = (aspect=1, title="Vorticity at t=10", xlabel="x", ylabel="y")
heatmap(ζ, colormap=:balance, colorrange=(-ζlim, ζlim); axis)

save("first_example.pdf", current_figure())
