using Oceananigans

"""
    advect_tracer(tracer_advection, closure=nothing, Nx=128)

Advect a tracer with a top-hat profile until `stop_time` using `Nx` grid points.
"""
function advect_tracer(tracer_advection; closure=nothing, Nx=128, stop_time=1)
    grid = RectilinearGrid(size=Nx, x=(-4, 8), halo=6, topology=(Periodic, Flat, Flat))

    velocities = PrescribedVelocityFields(u=1)
    model = HydrostaticFreeSurfaceModel(; grid, closure, velocities,
                                        tracer_advection, tracers = :c)

    set!(model, c= x -> abs(x) > 1 ? 0 : 1)
    simulation = Simulation(model; Δt=0.1/Nx, stop_time)
    run!(simulation)
    return model.tracers.c
end

using GLMakie

set_theme!(Theme(fontsize=18, linewidth=2, linealpha=0.6))
fig = Figure(size=(1000, 300))
ax = Axis(fig[1, 1], xlabel="x") #, ylabel="c(t=4)")

stop_time = 0.1
c = advect_tracer(Centered(order=2); stop_time)
lines!(ax, c, label="t=0.1, Centered(order=2)")

c = advect_tracer(UpwindBiased(order=3); stop_time)
lines!(ax, c, label="t=0.1, UpwindBiased(order=3)")

#c = advect_tracer(UpwindBiased(order=5); stop_time)
#lines!(ax, c, label="UpwindBiased(order=5)")

stop_time = 4

for N = (5, 11)
    scheme = WENO(order=N)
    c = advect_tracer(scheme; stop_time)
    lines!(ax, c, label="t=4, WENO(order=$N)")
end

solution(x) = abs(x) > 1 ? 0 : 1 
set!(c, solution)
lines!(ax, c, linestyle=:dash, color=(:black, 0.2), label="Initial condition")

solution(x) = abs(x-stop_time) > 1 ? 0 : 1 
set!(c, solution)
lines!(ax, c, linestyle=:dash, color=(:black, 0.6), label="Exact solution at t=4")

# c = advect_tracer(WENO(order=11), closure=ScalarDiffusivity(κ=0.04))
# lines!(ax, c, label="WENO(order=11) with κ=0.04")

xlims!(ax, -2, 6)
hidespines!(ax, :t, :r, :l)
hideydecorations!(ax, grid=false)

Legend(fig[1, 0], ax, framevisible=false)

display(fig)

save("advection_schemes.png", fig)

