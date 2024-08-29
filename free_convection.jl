using Oceananigans
using Oceananigans.Units
using Printf

grid = RectilinearGrid(GPU(),
                       topology = (Periodic, Periodic, Bounded),
                       size = (512, 512, 512),
                       x = (0, 128),
                       y = (0, 128),
                       z = (-64, 0))

top_buoyancy_bc = FluxBoundaryCondition(1e-7) # m² s⁻³
buoyancy_bcs = FieldBoundaryConditions(top=top_buoyancy_bc)

model = NonhydrostaticModel(; grid,
                            closure = AnisotropicMinimumDissipation(),
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (; b=buoyancy_bcs))

N² = 1e-6
Δz = minimum_zspacing(grid)
bᵢ(x, y, z) = N² * z + 1e-1 * N² * Δz * randn()
set!(model, b=bᵢ)

simulation = Simulation(model, Δt=10.0, stop_time=2hours)
conjure_time_step_wizard!(simulation, IterationInterval(10), cfl=0.2)

u, v, w = model.velocities
progress(sim) = @info @sprintf("Iter: %d, time: %s, max|w|: %.2e m s⁻¹",
                               iteration(sim), prettytime(sim), maximum(abs, w))
add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)

using CairoMakie

fig = Figure(size=(1200, 600)) 
axw = Axis(fig[1, 1])
axb = Axis(fig[1, 2])
heatmap!(axw, view(model.velocities.w, :, :, grid.Nz), colormap=:balance)
heatmap!(axb, view(model.tracers.b, :, 1, :))
