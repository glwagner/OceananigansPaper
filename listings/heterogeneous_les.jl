using Oceananigans
using Oceananigans.Units
using Printf

domain = :big
arch = GPU()
L = 2048
N = 512

Ny = 2N
Nz = 64
y = (-L, L)
z = (-256, 0)

if domain == :small
    Nx = Int(N/2)
    x = (0, L/2)
elseif domain == :big
    Nx = 4N
    x = (0, 4L)
end

f = 1e-4
τx = -1e-5
N² = 1e-6
J₀ = 1e-8
La = 0.3
Ri = 1.0
h₀ = 64
Δy = 128
σ = 2π / 7 # s
g = 9.81

prefix = "uniform_forcing_multi_scale_les_$(domain)_$N"

M² = sqrt(f^2 * N² / Ri)
Δb = M² * Δy
k = σ^2 / g
u★ = sqrt(abs(τx))
uˢ = u★ / La^2

@inline ramp(y, d=1) = max(y/d, zero(y))
@inline step(y, d=1) = max(-one(y), min(y / d, one(y)))

grid = RectilinearGrid(arch; size=(Nx, Ny, Nz), halo=(5, 5, 5),
                       x, y, z, topology=(Periodic, Bounded, Bounded))

@inline buoyancy_flux(x, y, t, p) = p.J₀ * (y > 0)
@inline wind_stress(x, y, t, p) = - p.τx * (y > 0)

#top_b_bc = FluxBoundaryCondition(buoyancy_flux, parameters=(; J₀))
top_b_bc = FluxBoundaryCondition(J₀)
top_u_bc = FluxBoundaryCondition(wind_stress, parameters=(; τx))

boundary_conditions = (; b = FieldBoundaryConditions(top=top_b_bc))
#                       u = FieldBoundaryConditions(top=top_u_bc))

coriolis = FPlane(; f)
advection = WENO(order=9)
closure = nothing

@inline ∂z_uˢ(z, t, p) = 2 * p.k * p.uˢ * exp(2 * p.k * z)
stokes_drift = UniformStokesDrift(; ∂z_uˢ, parameters=(; uˢ, k))

model = NonhydrostaticModel(; grid, advection, closure, coriolis,
                            boundary_conditions, # stokes_drift,
                            tracers = :b, buoyancy = BuoyancyTracer())

Δz = minimum_zspacing(grid)
db = N² * Δz
u★ = sqrt(abs(τx))
bᵢ(x, y, z) = N² * z - Δb / 2 * step(y, Δy) * ramp(z/h₀ + 1) + 1e-1 * db * (2rand() - 1)
uᵢ(x, y, z) = 1e-3 * u★ * (2rand() - 1)
set!(model, b=bᵢ, u=uᵢ, v=uᵢ, w=uᵢ)

Δt = 1minutes #0.1 * Δz / u★
simulation = Simulation(model; Δt, stop_time=30days)
conjure_time_step_wizard!(simulation, cfl=0.7, max_Δt=20minutes)

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

yzwriter = JLD2OutputWriter(model, outputs,
                            filename = prefix * "_yz.jld2",
                            indices = (Nx, :, :),
                            schedule = TimeInterval(30minutes),
                            overwrite_existing = true)

xzwriter = JLD2OutputWriter(model, outputs,
                            filename = prefix * "_xz.jld2",
                            indices = (:, 1, :),
                            schedule = TimeInterval(30minutes),
                            overwrite_existing = true)

w² = w^2
to_avg = (; u, v, w, w², b)
avgoutputs = NamedTuple(name => Average(to_avg[name], dims=1) for name in keys(to_avg))
avgwriter = JLD2OutputWriter(model, avgoutputs,
                             filename = prefix * "_averages.jld2",
                             schedule = TimeInterval(30minutes),
                             overwrite_existing = true)

checkpointer = Checkpointer(model; prefix,
                            schedule = TimeInterval(2days),
                            overwrite_existing = true,
                            cleanup = true)

simulation.output_writers[:xy] = xywriter
simulation.output_writers[:xz] = xzwriter
simulation.output_writers[:yz] = yzwriter
simulation.output_writers[:avg] = avgwriter
simulation.output_writers[:chk] = checkpointer

run!(simulation)

