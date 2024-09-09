using Oceananigans
using Oceananigans.BuoyancyModels
using SeawaterPolynomials.TEOS10
using Printf
using Statistics

grid = RectilinearGrid(topology = (Periodic, Flat, Bounded),
                       size = (128, 128), x = (0, 1), z = (-1, 0))

#equation_of_state = TEOS10.TEOS10EquationOfState()
equation_of_state = LinearEquationOfState(thermal_expansion=2e-4, haline_contraction=8e-5)
buoyancy = SeawaterBuoyancy(; equation_of_state)
#hydrostatic_pressure_anomaly = nothing
hydrostatic_pressure_anomaly = CenterField(grid)

model = NonhydrostaticModel(; grid, buoyancy,
                            advection = UpwindBiased(order=1),
                            hydrostatic_pressure_anomaly,
                            #timestepper = :RungeKutta3,
                            timestepper = :QuasiAdamsBashforth2,
                            tracers = (:T, :S))

ξᵢ(x, z) = 1e-9 * randn()
set!(model, T=20, S=35, u=ξᵢ, w=ξᵢ)

u, v, w = model.velocities
b = BuoyancyModels.BuoyancyField(model)
e = @at (Center, Center, Center) (u^2 + v^2 + w^2) / 2

simulation = Simulation(model, Δt=100, stop_iteration=4000)

B = Field(Average(b))
b′² = Field((b - B)^2)

function progress(sim)
    compute!(b′²)
    b²max = maximum(b′²)
    ebar = mean(e)
    msg = @sprintf("Iter: %d, max(b′²): %.2e, mean(e): %.2e",
                   iteration(sim), b²max, ebar)
    @info msg
    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))


outputs = merge(model.velocities, model.tracers, (; b, b′², e))
filename = "stable_motion_machine.jld2"

simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs; filename,
                                                    schedule = IterationInterval(10),
                                                    with_halos = true,
                                                    overwrite_existing = true)

run!(simulation)

ws = FieldTimeSeries(filename, "w")
bs = FieldTimeSeries(filename, "b")
b²s = FieldTimeSeries(filename, "b′²")
es = FieldTimeSeries(filename, "e")

Nt = length(es)
Es = zeros(Nt)
B²s = zeros(Nt)
for n = 1:Nt
    Es[n] = mean(es[n])
    B²s[n] = mean(b²s[n])
end

using GLMakie

fig = Figure(size=(600, 800))
axw = Axis(fig[3, 1], aspect=1, xlabel="x", ylabel="z", title="Vertical velocity")
axb = Axis(fig[3, 2], aspect=1, ylabel="z", xlabel="Buoyancy")
axe = Axis(fig[1, 1:2], ylabel="Mean kinetic energy", yscale=log10, xlabel="Time")
axb² = Axis(fig[2, 1:2], ylabel="Mean buoyancy variance", yscale=log10, xlabel="Time")

lines!(axe, es.times, Es)
lines!(axb², es.times, B²s)

slider = Slider(fig[4, 1:2], startvalue=Nt, range=1:Nt)
n = slider.value

tn = @lift es.times[$n]
En = @lift Es[$n]
B²n = @lift B²s[$n]

scatter!(axe, tn, En, marker=:circle, color=(:red, 0.5), markersize=20)
scatter!(axb², tn, B²n, marker=:circle, color=(:red, 0.5), markersize=20)

wn = @lift ws[$n]
bn = @lift bs[$n]

bmax = maximum(bs[end])
bmin = minimum(bs[end])
blims = (bmin, bmax)

heatmap!(axw, wn, colormap=:balance)
heatmap!(axb, bn) #, colorrange=blims, colormap=:balance)

display(fig)

record(fig, "stable_motion_machine.mp4", 1:Nt, framerate=24) do nn
    n[] = nn
end

#=
#ρ = seawater_density(model)
#ρs = FieldTimeSeries(filename, "ρ")
#axρ = Axis(fig[2, 3], aspect=1, ylabel="z", xlabel="Horizontally-averaged density")
#ρn = @lift ρs[$n]
ρt = CenterField(ρs.grid)
R = Field(Average(ρt, dims=1))
Rn = @lift begin
    parent(ρt) .= parent(ρs[$n])
    compute!(R)
    interior(R, 1, 1, :)
end
#lines!(axρ, Rn, z)
=#


