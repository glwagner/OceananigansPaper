using Oceananigans
using Oceananigans.Units
using Oceananigans.BuoyancyModels

grid = RectilinearGrid(topology = (Periodic, Flat, Bounded),
                       size = (128, 128), x = (0, 100), z = (-100, 0))

equation_of_state = LinearEquationOfState(thermal_expansion=2e-4, haline_contraction=8e-5)
buoyancy = SeawaterBuoyancy(; equation_of_state)
hydrostatic_pressure_anomaly = nothing
#hydrostatic_pressure_anomaly = CenterField(grid)

model = NonhydrostaticModel(; grid, buoyancy, advection = WENO(),
                            hydrostatic_pressure_anomaly,
                            timestepper = :RungeKutta3, tracers = (:T, :S))

Tᵢ(x, z) = 20 + 1e-5 * randn()
Sᵢ(x, z) = 35 + 1e-5 * randn()
set!(model, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=10minutes, stop_time=100days)
progress(sim) = @info string("Iter: ", iteration(sim), ", max|w|: ", maximum(abs, w))
add_callback!(simulation, progress, IterationInterval(10))

u, v, w = model.velocities
b = BuoyancyModels.buoyancy(model)
e = @at (Center, Center, Center) (u^2 + v^2 + w^2) / 2

outputs = merge(model.velocities, model.tracers, (; b, e))
filename = "test_linear_EOS_better_perturbation.jld2"

simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs; filename,
                                                    schedule = IterationInterval(10),
                                                    with_halos = true,
                                                    overwrite_existing = true)

run!(simulation)

ws = FieldTimeSeries(filename, "w")
bs = FieldTimeSeries(filename, "b")
es = FieldTimeSeries(filename, "e")

Nt = length(es)
Es = zeros(Nt)
for n = 1:Nt
    Es[n] = sum(es[n])
end

using GLMakie

fig = Figure(size=(1200, 600))
axw = Axis(fig[2, 1], aspect=1, xlabel="x", ylabel="z", title="Vertical velocity")
axb = Axis(fig[2, 2], aspect=1, ylabel="z", xlabel="Horizontally-averaged buoyancy")
axe = Axis(fig[1, 1:3], ylabel="Kinetic energy", xlabel="Time")

lines!(axe, es.times, Es)

slider = Slider(fig[3, 1:3], startvalue=Nt, range=1:Nt)
n = slider.value

wn = @lift ws[$n]
bn = @lift bs[$n]

heatmap!(axw, wn)

bt = CenterField(bs.grid)
B = Field(Average(bt, dims=1))

Bn = @lift begin
    parent(bt) .= parent(bs[$n])
    compute!(B)
    interior(B, 1, 1, :)
end

z = znodes(bs)
lines!(axb, Bn, z)

display(fig)

