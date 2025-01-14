using Oceananigans
using Oceananigans.Units
using Printf

arch = GPU()
x = y = (0, 4096)
z = (-256, 0)
Nx = Ny = 512
Nz = 64

grid = RectilinearGrid(arch, size=(Nx, Ny, Nz); x, y, z, topology=(Periodic, Periodic, Bounded))

# Background flow is Ψ = Λ y z
# corresponding to geostrophic shear U = Λ z
# and buoyancy B = - f Λ y satisfying thermal wind.
# Ri = N² / u_z² = N² / Λ²
# so Λ = N / √(Ri)

f = 1e-4
N² = 1e-6
Ri = 1
Λ = sqrt(N² / Ri)
parameters = (; Λ, f)

@inline uᵇᵍ(x, y, z, t, p) = + p.Λ * z
@inline bᵇᵍ(x, y, z, t, p) = - p.f * p.Λ * y
background_fields = (u = BackgroundField(uᵇᵍ; parameters),
                     b = BackgroundField(bᵇᵍ; parameters))

model = NonhydrostaticModel(; grid, background_fields,
                            coriolis = FPlane(; f),
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            advection = WENO(order=5))

Δz = minimum_zspacing(grid)
bᵢ(x, y, z) = N² * z + 1e-2 * N² * Δz * (2rand() - 1)
set!(model, b=bᵢ)

simulation = Simulation(model, Δt=10minutes, stop_time=10days)
conjure_time_step_wizard!(simulation, cfl=0.7)

wallclock = Ref(time_ns())
function progress(sim)
    elapsed = 1e-9 * (time_ns() - wallclock[])
    u, v, w = sim.model.velocities
    msg = @sprintf("Iter: %d, time: %s, Δt: %s, wall time: %s, max|𝐮|: (%.2f, %.2f, %.2f)",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(elapsed),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w)) 
    @info msg
    wallclock[] = time_ns()
    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)
s = @at (Center, Center, Center) sqrt(u^2 + v^2)
∇b² = @at (Center, Center, Center) ∂x(b)^2 + ∂y(b)^2
outputs = merge(model.velocities, model.tracers, (; ζ, s, ∇b²))

xy = JLD2OutputWriter(model, outputs,
                      filename = "eady_les.jld2",
                      schedule = TimeInterval(10minutes),
                      indices = (:, :, Nz),
                      overwrite_existing = true)

simulation.output_writers[:xy] = xy

run!(simulation)

