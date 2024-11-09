using Oceananigans
using Oceananigans.Units
using Printf

x = y = (0, 512)
z = (-256, 0)
grid = RectilinearGrid(CPU(); size=(32, 32, 32), halo=(5, 5, 5), x, y, z)

@inline τx(x, y, t, p) = p.τ₀ * exp(-t^2 / (2 * p.T^2))
u_top_bc = FluxBoundaryCondition(τx, parameters=(τ₀=1e-3, T=12hours))
u_bcs = FieldBoundaryConditions(top=u_top_bc)

model = NonhydrostaticModel(; grid,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (; u=u_bcs),
                            advection = WENO(order=9),
                            coriolis = FPlane(f=1e-4))
                            
N² = 1e-5
Δz = minimum_zspacing(grid)
bᵢ(x, y, z) = 1e-5 * z + 1e-2 * N² * Δz * randn()
set!(model, b=bᵢ)

simulation = Simulation(model, Δt=1minute, stop_time=24hours)
conjure_time_step_wizard!(simulation, cfl=0.7)

function progress(sim)
    u, v, w = sim.model.velocities
    msg = @sprintf("Iter: %d, time: %s, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w))
    @info msg
    return nothing
end
                   
add_callback!(simulation, progress, IterationInterval(10))                   

Nx, Ny, Nz = size(grid)
outputs = merge(model.velocities, model.tracers)
filename = "pulse_of_wind"

xyow = JLD2OutputWriter(model, outputs; filename=filename * "_xy.jld2",
                        schedule=TimeInterval(10minutes), indices=(:, :, Nz), overwrite_existing=true)

yzow = JLD2OutputWriter(model, outputs; filename=filename * "_yz.jld2",
                        schedule=TimeInterval(10minutes), indices=(1, :, :), overwrite_existing=true)

xzow = JLD2OutputWriter(model, outputs; filename=filename * "_xz.jld2",
                        schedule=TimeInterval(10minutes), indices=(:, 1, :), overwrite_existing=true)
                        
simulation.output_writers[:xy] = xyow
simulation.output_writers[:xz] = xzow
simulation.output_writers[:yz] = yzow

run!(simulation)

