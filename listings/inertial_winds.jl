using Oceananigans
using Oceananigans.Units

grid = RectilinearGrid(GPU();
                       size = (512, 512, 256),
                       halo = (5, 5, 5),
                       x = (0, 512),
                       y = (0, 512),
                       z = (-256, 0))

@inline τx(x, y, t, δ) = - 1e-3 * exp(-(t - 2δ)^2 / 2δ^2)
u_top_bc = FluxBoundaryCondition(τx, parameters=6hours)
u_bcs = FieldBoundaryConditions(top=u_top_bc)

@inline ∂z_uˢ(z, t) = 0.01 * exp(z / 8)
stokes_drift = UniformStokesDrift(; ∂z_uˢ)

@inline mask(z, w, H) = max(zero(z), 1 - (z + H) / w)
@inline sponge(x, y, z, t, q, ω) = - ω * q * mask(z, 32, 256)
u_forcing = Forcing(sponge, field_dependencies=:u, parameters=0.05)
v_forcing = Forcing(sponge, field_dependencies=:v, parameters=0.05)
w_forcing = Forcing(sponge, field_dependencies=:w, parameters=0.05)

model = NonhydrostaticModel(; grid, stokes_drift,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            forcing = (u=u_forcing, v=v_forcing, w=w_forcing),
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

wall_time = Ref(time_ns())
sim_time = Ref(time(simulation))

function progress(sim)
    u, v, w = sim.model.velocities
    elapsed = 1e-9 * (time_ns() - wall_time[])
    simulated = time(sim) - sim_time[]
    SDPD = simulated / elapsed

    msg = @sprintf("Iter: %d, time: %s, Δt: %s, elapsed: %s, SDPD: %.1f, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(elapsed), SDPD,
                   maximum(abs, u), maximum(abs, v), maximum(abs, w))
    @info msg
    wall_time[] = time_ns()
    sim_time[] = time(sim)
    return nothing
end
                   
add_callback!(simulation, progress, IterationInterval(10))                   

Nx, Ny, Nz = size(grid)
outputs = merge(model.velocities, model.tracers)
filename = string("wind_pulse_", grid.Nx)
save_interval = 5minutes

xyow = JLD2Writer(model, outputs; filename=filename * "_xy.jld2",
                  schedule=TimeInterval(save_interval), indices=(:, :, Nz), overwrite_existing=true)

yzow = JLD2Writer(model, outputs; filename=filename * "_yz.jld2",
                  schedule=TimeInterval(save_interval), indices=(1, :, :), overwrite_existing=true)

xzow = JLD2Writer(model, outputs; filename=filename * "_xz.jld2",
                  schedule=TimeInterval(save_interval), indices=(:, 1, :), overwrite_existing=true)
                        
simulation.output_writers[:xy] = xyow
simulation.output_writers[:xz] = xzow
simulation.output_writers[:yz] = yzow

run!(simulation)

