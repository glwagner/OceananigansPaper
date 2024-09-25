using Oceananigans

grid = RectilinearGrid(size=64, x=(-4, 8), topology=(Periodic, Flat, Flat)) 

model = NonhydrostaticModel(grid=grid, tracers=:c, advection=WENO(order=5))
cᵢ(x) = exp(-x^2 / 2)
set!(model, u=1, c=cᵢ)

using CairoMakie
lines(model.tracers.c, label="Tracer at t=0")

simulation = Simulation(model, Δt=0.01, stop_time=1)
run!(simulation)
lines!(model.tracers.c, label="Tracer at t=1")

