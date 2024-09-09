using Oceananigans

grid = RectilinearGrid(size = (256, 256),
                       x = (-π, π), y = (-π, π),
                       topology = (Periodic, Bounded, Flat))
                       
coriolis = BetaPlane(f₀=0, β=20)
model = NonhydrostaticModel(; grid, coriolis, advection=WENO())

ϵ(x, y) = 2rand() - 1
set!(model, u=ϵ, v=ϵ)

simulation = Simulation(model; Δt=0.01, stop_time=40)

using Printf

function progress(sim)
    msg = @sprintf("Iter: %d, time: %.2f", iteration(sim), time(sim))
    @info msg
    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)

u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u))
compute!(ζ)

U = Field(Average(u, dims=1))
compute!(U)

using GLMakie

fig = Figure(size=(800, 600))
axζ = Axis(fig[1, 1], aspect=1, xlabel="x", ylabel="x")
axU = Axis(fig[1, 2], xlabel="U", ylabel="y", yaxisposition=:right)

heatmap!(axζ, ζ, colormap=:balance, colorrange=(-ζlim, ζlim))
lines!(axU, interior(U, 1, :, 1), ynodes(U))
ylims!(axU, -π, π)

colsize!(fig.layout, 2, Relative(0.2))

save("beta_plane_turbulence.png", current_figure())

display(current_figure())
