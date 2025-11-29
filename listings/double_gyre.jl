using Oceananigans
using Oceananigans: WallTimeInterval
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using Oceananigans.Grids: ExponentialDiscretization
using Printf
using CUDA

arch = GPU()
Nx = 12 #* 22   FJP change these back
Ny = 12 #* 20
Nz = 10 #0
Lz = 4000
z_faces = ExponentialDiscretization(Nz, -Lz, 0; scale = Lz/5, bias = :right)
stop_time = 360days * 4

grid = LatitudeLongitudeGrid(arch,
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             longitude = (-11, 11),
                             latitude = (30, 50),
                             z = z_faces)

@show grid

δλ, δφ, H∞ = 10, 20, 6000
paraboloid(λ, φ) = λ^2 / δλ^2 + (φ - 30)^2 / δφ^2
bottom_height(λ, φ) = - H∞ * (1.1 - paraboloid(λ, φ))
grid = ImmersedBoundaryGrid(grid, PartialCellBottom(bottom_height))

@inline τx(λ, φ, t, τ₀) = -τ₀ * (1 - cos(2π * (φ - 30) / 20))
u_top_bc = FluxBoundaryCondition(τx, parameters=2e-4)
u_bcs = FieldBoundaryConditions(top=u_top_bc)

@inline b★(φ, p) = - p.Δb * sin(π/2 * (φ - 30) / 20)
@inline b_restoring(λ, φ, t, b, p) = p.ω * (b - b★(φ, p))
b_top_bc = FluxBoundaryCondition(b_restoring, field_dependencies=:b, parameters=(ω=1/1days, Δb=0.2))
b_bcs = FieldBoundaryConditions(top=b_top_bc)

model = HydrostaticFreeSurfaceModel(; grid,
                                    coriolis = HydrostaticSphericalCoriolis(),
                                    momentum_advection = WENOVectorInvariant(),
                                    tracer_advection = WENO(order=9),
                                    closure = CATKEVerticalDiffusivity(),
                                    tracers = (:b, :e),
                                    buoyancy = BuoyancyTracer(),
                                    boundary_conditions = (u=u_bcs, b=b_bcs))

N² = 2e-5
bᵢ(x, y, z) = N² * z
set!(model, b=bᵢ)

simulation = Simulation(model; Δt=1minutes, stop_time)
conjure_time_step_wizard!(simulation, cfl=0.2)

function progress(sim)
    u, v, w = sim.model.velocities
    @info @sprintf("Iter: %d, time: %s, Δt: %s, max|w|: %.2e",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), maximum(w))

    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

e = model.tracers.e
b = model.tracers.b
u, v, w = model.velocities
κc = model.closure_fields.κc
outputs = (; u, v, w, b, e, κc)

Nx, Ny, Nz = size(grid)
i = floor(Int, Nx/2)
j = floor(Int, Ny/2)
indiceses = [(:, 1, :), (i, :, :),
             (:, :, Nz), (:, :, 1), (:, :, Nz-8)]
names = [:xz, :yz, :xy1, :xy2, :xy3]
prefix = "double_gyre_Nx$(Nx)_Nz$(Nz)"

for (name, indices) in zip(names, indiceses)
    output_writer = JLD2Writer(model, outputs; indices,
                               filename = string(prefix, "_", name),
                               schedule = TimeInterval(2day),
                               with_halos = true,
                               overwrite_existing = true)
                                          
    simulation.output_writers[name] = output_writer
end

checkpointer = Checkpointer(model,
                            prefix = "double_gyre",
                            schedule = WallTimeInterval(3hours),
                            cleanup = true,
                            overwrite_existing = true)

simulation.output_writers[:chk] = checkpointer

run!(simulation)

