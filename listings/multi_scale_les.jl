using Oceananigans
using Oceananigans.Units
using Printf

arch = GPU()
Nx = 2048
Ny = 1024
Nz = 48
L = 2048
x = (-2L, 2L)
y = (-L, L)
z = (-192, 0)

f = 1e-4
τx = -1e-6
N² = 5e-6
J₀ = 5e-9
La = 0.3
Ri = 1.0
h₀ = 64
Δy = 128
σ = 2π / 7 # s
g = 9.81

prefix = "multi_scale_les_half_size_$Nx"

M² = sqrt(f^2 * N² / Ri)
Δb = M² * Δy
k = σ^2 / g
u★ = sqrt(abs(τx))
uˢ = u★ / La^2

@inline ramp(y, d=1) = max(y/d, zero(y))
@inline step(y, d=1) = max(-one(y), min(y / d, one(y)))

grid = RectilinearGrid(arch; size=(Nx, Ny, Nz), halo=(5, 5, 5),
                       x, y, z, topology=(Periodic, Bounded, Bounded))

@inline buoyancy_flux(x, y, t, p) = p.J₀ * (1 + step(y, p.Δy))

top_b_bc = FluxBoundaryCondition(buoyancy_flux, parameters=(; J₀, Δy))
top_u_bc = FluxBoundaryCondition(τx)

boundary_conditions = (b = FieldBoundaryConditions(top=top_b_bc),
                       u = FieldBoundaryConditions(top=top_u_bc))

coriolis = FPlane(; f)
advection = WENO(order=9)
closure = nothing

@inline ∂z_uˢ(z, t, p) = 2 * p.k * p.uˢ * exp(2 * p.k * z)
stokes_drift = UniformStokesDrift(; ∂z_uˢ, parameters=(; uˢ, k))

model = NonhydrostaticModel(; grid, advection, closure, coriolis,
                            boundary_conditions, stokes_drift,
                            tracers = :b, buoyancy = BuoyancyTracer())

Δz = minimum_zspacing(grid)
db = N² * Δz
u★ = sqrt(abs(τx))
bᵢ(x, y, z) = N² * z - Δb / 2 * step(y, Δy) * ramp(z/h₀ + 1) + 1e-1 * db * (2rand() - 1)
uᵢ(x, y, z) = 1e-3 * u★ * (2rand() - 1)
set!(model, b=bᵢ, u=uᵢ, v=uᵢ, w=uᵢ)

Δt = 1minutes #0.1 * Δz / u★
simulation = Simulation(model; Δt, stop_time=30days)
conjure_time_step_wizard!(simulation, cfl=0.7, max_Δt=10minutes)

stopwatch = Ref(time_ns())
function progress(sim)
    elapsed = 1e-9 * (time_ns() - stopwatch[])
    msg = @sprintf("Iter: %d, time: %s, Δt: %s, wall time: %s",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(elapsed))

    u, v, w = sim.model.velocities
    msg *= @sprintf(", max|u|: (%.2e, %.2e, %.2e)",
                    maximum(abs, interior(u)),
                    maximum(abs, interior(v)),
                    maximum(abs, interior(w)))

    @info msg

    stopwatch[] = time_ns()

    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

b = model.tracers.b
u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)
outputs = (; u, v, w, ζ, b)
xywriter = JLD2OutputWriter(model, outputs,
                            filename = prefix * "_xy.jld2",
                            indices = (:, :, Nz),
                            schedule = TimeInterval(30minutes),
                            overwrite_existing = true)

w² = w^2
to_avg = (; u, v, w, w², b)
avgoutputs = NamedTuple(name => Average(to_avg[name], dims=1) for name in keys(to_avg))
avgwriter = JLD2OutputWriter(model, avgoutputs,
                             filename = prefix * "_averages.jld2",
                             schedule = TimeInterval(30minutes),
                             overwrite_existing = true)

checkpointer = Checkpointer(model,
                            prefix = "multi_scale_les",
                            schedule = TimeInterval(14days),
                            cleanup = true)

simulation.output_writers[:xy] = xywriter
simulation.output_writers[:avg] = avgwriter

run!(simulation)

