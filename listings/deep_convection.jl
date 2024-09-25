using Oceananigans
using Oceananigans.Units
using Printf

Ω = 7.292115e-5
φ = 50
fz = 2Ω * cos(φ)
fy = 2Ω * sin(φ)
L = 512

grid = RectilinearGrid(GPU(),
                       topology = (Periodic, Periodic, Bounded),
                       size = (384, 384, 384),
                       halo = (6, 6, 6),
                       x = (-L, L),
                       y = (-L, L),
                       z = (-2L, 0))

@inline circle(x, y) = (x^2 + y^2) < 200^2
@inline Jb(x, y, t) = 1e-6 * circle(x, y)
top_b_bc = FluxBoundaryCondition(Jb)
b_bcs = FieldBoundaryConditions(top=top_b_bc)

u★ = 3e-3
top_u_bc = FluxBoundaryCondition(-u★^2)
u_bcs = FieldBoundaryConditions(top=top_u_bc)

La = 0.3
@inline ∂z_uˢ(z, t, (k, Uˢ)) = Uˢ / 2k * exp(2k * z)
stokes_drift = UniformStokesDrift(; ∂z_uˢ, parameters = (2π/100, u★/La^2))

model = NonhydrostaticModel(; grid, # stokes_drift,
                            advection = WENO(order=9),
                            tracers = :b,
                            coriolis = ConstantCartesianCoriolis(; fx=0, fy, fz),
                            timestepper = :RungeKutta3,
                            buoyancy = BuoyancyTracer(),
                            boundary_conditions = (b=b_bcs, u=u_bcs))

N² = 1e-6
Δz = minimum_zspacing(grid)
bᵢ(x, y, z) = N² * z
wᵢ(x, y, z) = 1e-3 * randn()
set!(model, b=bᵢ, w=wᵢ)

simulation = Simulation(model, Δt=1.0, stop_time=5*72hours)
conjure_time_step_wizard!(simulation, IterationInterval(3), cfl=0.7, max_Δt=20)

u, v, w = model.velocities
progress(sim) = @info @sprintf("Iter: %d, time: %s, Δt: %s, max|u|: %.2e m s⁻¹, max|w|: %.2e m s⁻¹",
                               iteration(sim), prettytime(sim), prettytime(sim.Δt),
                               maximum(abs, interior(u)),
                               maximum(abs, interior(w)))

add_callback!(simulation, progress, IterationInterval(10))

b = model.tracers.b
outputs = (; u, v, w, b)
Nx, Ny, Nz = size(grid)
i = floor(Int, Nx/2)
j = floor(Int, Ny/2)
indiceses = [(:, j, :), (:, :, Nz), (i, :, :), (:, 1, :), (1, :, :)]
names = [:xz, :xy, :yz, :xz1, :yz1]
prefix = "deep_convection_no_stokes_$Nz"

for (name, indices) in zip(names, indiceses)
    output_writer = JLD2OutputWriter(model, outputs; indices,
                                     filename = string(prefix, "_", name),
                                     schedule = TimeInterval(2minutes),
                                     with_halos = true,
                                     overwrite_existing = true)
                                          
    simulation.output_writers[name] = output_writer
end

checkpointer = Checkpointer(model,
                            prefix = "checkpoint_$prefix",
                            schedule = TimeInterval(24hours),
                            overwrite_existing = true,
                            cleanup = true)
                                          
@show model
@show simulation

run!(simulation)

