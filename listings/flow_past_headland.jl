using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials

H, L, δ, N = 256, 1024, 512, 64
x, y, z = (-2L, 2L), (-L, L), (-H, 0)

grid = RectilinearGrid(GPU(); size=(2N, 4N, N), halo=(6, 6, 6),
                       x, y, z, topology=(Periodic, Bounded, Bounded))

bowl(y) = 1 + (y / L)^2
wedge(x, y) = 1 + (y + abs(x)) / δ
bowl_wedge(x, y) = -H * max(bowl(y), wedge(x, y))
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bowl_wedge))

T₂ = 12.421hours
U₂ = 0.15 # m/s
@inline Fu(x, y, z, t, p) = 2π * p.U₂ / p.T₂ * cos(2π * t / p.T₂)

equation_of_state = SeawaterPolynomials.TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

model = NonhydrostaticModel(; grid, buoyancy, advection = WENO(order=9),
                            forcing = (; u=Forcing(Fu, parameters=(; U₂, T₂))),
                            coriolis = FPlane(latitude=47.5),
                            tracers = (:T, :S))

Tᵢ(x, y, z) = 12 + 4z / H
set!(model, T=Tᵢ, S=32, u=0.15)

simulation = Simulation(model, Δt=10, stop_time=4days)
conjure_time_step_wizard!(simulation, cfl=0.7)

using Printf

wallclock = Ref(time_ns())

function progress(sim)
    u, v, w = sim.model.velocities
    ΔT = 1e-9 * (time_ns() - wallclock[])
    msg = @sprintf("(%d) t: %s, Δt: %s, wall time: %s, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(ΔT),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w))
    @info msg
    wallclock[] = time_ns()
    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)
s = @at (Center, Center, Center) sqrt(u^2 + v^2)
xy = JLD2OutputWriter(model, merge(model.velocities, model.tracers, (; ζ, s)),
                      indices = (:, :, grid.Nz),
                      schedule = TimeInterval(10minutes),
                      filename = "flow_past_headland.jld2",
                      overwrite_existing = true)

simulation.output_writers[:xy] = xy

run!(simulation)

