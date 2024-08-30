#=
using Oceananigans
using Oceananigans.Models: seawater_density
using SeawaterPolynomials: TEOS10EquationOfState

grid = RectilinearGrid(GPU(),
                       size = (2048, 2048),
                       x = (-0.25, 0.25),
                       z = (-0.5, 0.0),
                       topology = (Bounded, Flat, Bounded))

closure = ScalarDiffusivity(ν=1.15e-6, κ=(T=1e-7, S=1e-9))

equation_of_state = TEOS10EquationOfState(reference_density=1000)
buoyancy = SeawaterBuoyancy(; equation_of_state, gravitational_acceleration=9.81) 

model = NonhydrostaticModel(; grid, buoyancy, closure, tracers=(:T, :S), timestepper=:RungeKutta3)
                            
T₁, T₂ = 1, 7.55 # ᵒC
Tᵢ(x, z) = z <= -0.25 ? T₁ : T₂
Ξᵢ(x, z) = 1e-2 * randn()
set!(model, T=Tᵢ, S=0, u=Ξᵢ, v=Ξᵢ, w=Ξᵢ) 

Δx = minimum_xspacing(grid)
Δt = 0.2 * Δx^2 / closure.ν
simulation = Simulation(model; Δt, stop_time=1000)

u, v, w = model.velocities
progress(sim) = @info string("Iter: ", iteration(sim),
                             ", time: ", prettytime(sim),
                             ", max(w): ", maximum(abs, w))

add_callback!(simulation, progress, IterationInterval(10))

ρ = Field(seawater_density(model))
T = model.tracers.T

output_writer = JLD2OutputWriter(model, (; ρ, T),
                                 filename = "cabbeling",
                                 schedule = TimeInterval(2),
                                 overwrite_existing = true)
                                        
simulation.output_writers[:jld2] = output_writer

run!(simulation)
=#

using CairoMakie

ρt = FieldTimeSeries("cabbeling.jld2", "ρ")
Tt = FieldTimeSeries("cabbeling.jld2", "T")

Nt = length(ρt)
n = Observable(1)
ρ = @lift ρt[$n]
T = @lift Tt[$n]

ρ = Field(seawater_density(model))
compute!(ρ)

fig = Figure(size=(700, 400))

axρ = Axis(fig[2, 1], aspect=1, xlabel="x (m)", ylabel="z (m)")
hm = heatmap!(axρ, ρ, colormap=CairoMakie.Reverse(:grays), colorrange=(999.90, 999.98))
Colorbar(fig[1, 1], hm, label="Density (kg m⁻³)", vertical=false)

axT = Axis(fig[2, 2], aspect=1, xlabel="x (m)", ylabel="z (m)")
hm = heatmap!(axT, T, colormap=:thermal, colorrange=(1.2, 7.3))
Colorbar(fig[1, 2], hm, label="Temperature (ᵒC)", vertical=false)

record(fig, "cabbeling.mp4", 1:Nt, framerate=12) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

#=
#save("cabbeling.png", current_figure())
ρmax = maximum(ρ)
ρmin = minimum(ρ)
ρmid = (ρmax + ρmin) / 2
Δρ = ρmax - ρmin
ρhi = ρmid + 0.4 * Δρ
ρlo = ρmid - 0.4 * Δρ

using CairoMakie
fig = Figure(size=(600, 1200))
axT = Axis(fig[1, 1])
axρ = Axis(fig[1, 2])
heatmap!(axT, model.tracers.T, colormap=:thermal)
heatmap!(axρ, ρ, colormap=:grays, colorrange=(ρlo, ρhi))

save("cabbeling.png", current_figure())
=#
