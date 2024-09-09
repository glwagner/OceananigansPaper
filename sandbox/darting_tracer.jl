using Oceananigans

# A two-dimensional square grid on the CPU.
grid = RectilinearGrid(CPU(),
                       size = (256, 256),
                       x = (-π, π),
                       y = (-π, π),
                       topology = (Periodic, Periodic, Flat))

# Idea: make the tracer move from one point to another
# Update the new point randomly with a callback

# Updatable positions
X₀ = Ref((0.0, 0.0))
X₁ = Ref((0.0, 0.0))

function darting_source(x, y, t)
    ϵ = 1e-6
    x₀, y₀ = X₀[]
    x₁, y₁ = X₁[]
    
    # Compute new position
    x₀ = (1 - ϵ) * x₀ + ϵ * x₁
    y₀ = (1 - ϵ) * y₀ + ϵ * y₁
    X₀[] = (x₀, y₀)
    
    # Compute source term
    δ = 0.1
    return ((x - x₀)^2 + (y - y₀)^2) < δ^2
end
                       
# Build a model with the "Weighted Essentially Non-Oscillatory" (WENO)
# advection scheme, which introduces implicit numerical diffusion.
model = NonhydrostaticModel(; grid, advection=WENO(),
                            tracers=:c, forcing=(; c=darting_source))

# Set a random initial condition and then run the simulation. 
ϵ(x, y) = 2rand() - 1
set!(model, u=ϵ, v=ϵ)

simulation = Simulation(model; Δt=0.001, stop_time=1)

using Printf

function print_progress(sim)
    u, v, w = sim.model.velocities
    max_u = max(maximum(abs, u), maximum(abs, v))
    @info @sprintf("Iteration: %d, time: %.1f, max|u|: %.2e, (x₀, y₀): (%.2f, %.2f)",
                   iteration(sim), time(sim), max_u, X₀[][1], X₀[][2]) 
    return nothing
end

add_callback!(simulation, print_progress, IterationInterval(10))

function update_target(sim)
    x₁ = 2π * (rand() - 0.5)
    y₁ = 2π * (rand() - 0.5)
    X₁[] = (x₁, y₁)
    return nothing
end

add_callback!(simulation, update_target, TimeInterval(0.1))

run!(simulation)

using CairoMakie

clim = maximum(model.tracers.c) / 2
heatmap(model.tracers.c, colormap=:thermal, colorrange=(0, clim), axis=(; aspect=1))

display(current_figure())
save("second_example.pdf", current_figure())
