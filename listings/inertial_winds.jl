using Oceananigans
using Oceananigans.Units

grid = RectilinearGrid(GPU();
                       size = (256, 256, 256),
                       halo = (5, 5, 5),
                       x = (0, 512),
                       y = (0, 512),
                       z = (-256, 0))

@inline τx(x, y, t, δ) = - 1e-3 * exp(-t^2 / 2δ)
u_top_bc = FluxBoundaryCondition(τx, parameters=12hours)
u_bcs = FieldBoundaryConditions(top=u_top_bc)

@inline ∂z_uˢ(z, t) = 0.01 * exp(z / 8)
stokes_drift = UniformStokesDrift(; ∂z_uˢ)

model = NonhydrostaticModel(; grid, stokes_drift,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (; u=u_bcs),
                            advection = WENO(order=9),
                            coriolis = FPlane(f=1e-4))
                            
N² = 1e-5
Δz = minimum_zspacing(grid)
bᵢ(x, y, z) = N² * z + 1e-2 * N² * Δz * randn()
set!(model, b=bᵢ)


simulation = Simulation(model, Δt=10.0, stop_time=36hours)
conjure_time_step_wizard!(simulation, cfl=0.7)

using Printf
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

