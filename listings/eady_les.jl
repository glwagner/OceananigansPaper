using Oceananigans
using Oceananigans.Units
using Printf

arch = CPU()
x = y = (0, 4kilometers)
z = (-256, 0)
Nx = Ny = 64
Nz = 16

grid = RectilinearGrid(arch, size=(Nx, Ny, Nz); x, y, z)

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

simulation = Simulation(model, Δt=10minutes, stop_time=30days)
conjure_time_step_wizard!(simulation, cfl=0.7)

wallclock = Ref(time_ns())
function progress(sim)
    elapsed = 1e-9 * (time_ns() - wallclock[])
    msg = @sprintf("Iter: %d, time: %s, Δt: %s, wall time: %s",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(elapsed))
    @info msg
    wallclock[] = time_ns()
    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))


run!(simulation)

using GLMakie

u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u), indices=(:, :, Nz))
compute!(ζ)

heatmap(ζ)
