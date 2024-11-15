using Oceananigans
using Oceananigans.Units

L = 512
H = 128
δ = L/2
N = 8

x = (-2L, 2L)
y = (-L, L)
z = (-H, 0)

grid = RectilinearGrid(CPU(); size=(8N, 16N, 2N), halo=(6, 6, 6),
                       x, y, z, topology=(Periodic, Bounded, Bounded))

bowl(y) = -H * (1 + (y/L)^2)
wedge(x, y) = y > -L + δ - abs(x)
bowl_wedge(x, y) = bowl(y) * wedge(x, y)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bowl_wedge))

model = NonhydrostaticModel(; grid,
                            advection = WENO(order=5),
                            coriolis = FPlane(f=1e-4),
                            tracers = :b,
                            buoyancy = BuoyancyTracer())

bᵢ(x, y, z) = 1e-5 * z
set!(model, b=bᵢ, u=1)

simulation = Simulation(model, Δt=10, stop_iteration=100)
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

using GLMakie

ut = FieldTimeSeries("flow_past_headland.jld2", "u")
Nt = length(ut)

fig = Figure(size=(1200, 600))
axu = Axis(fig[1, 1], aspect=2)
slider = Slider(fig[2, 1], range=1:Nt, startvalue=1)
n = slider.value

un = @lift ut[$n]
heatmap!(axu, un)

display(fig)

