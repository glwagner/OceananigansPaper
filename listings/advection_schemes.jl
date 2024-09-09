using Oceananigans

"""
    advect_tracer(tracer_advection, closure=nothing, Nx=128)

Advect a tracer with a top-hat profile until `stop_time` using `Nx` grid points.
"""
function advect_tracer(tracer_advection; closure=nothing, Nx=128)
    grid = RectilinearGrid(size=Nx, x=(-4, 8), halo=7, topology=(Periodic, Flat, Flat))

    velocities = PrescribedVelocityFields(u=1)
    model = HydrostaticFreeSurfaceModel(; grid, closure, velocities,
                                        tracer_advection, tracers = :c)

    set!(model, c= x -> abs(x) > 1 ? 0 : 1)
    simulation = Simulation(model; Δt=0.1 * 12/Nx, stop_time=4)
    run!(simulation)
    return model.tracers.c
end

using GLMakie
fig = Figure(size=(600, 300))
ax = Axis(fig[1, 1], xlabel="x", ylabel="c(t=4)")

for N = (5, 11)
    scheme = WENO(order=N)
    c = advect_tracer(scheme)
    lines!(ax, c, label="WENO(order=$N)")
end

c = advect_tracer(UpwindBiased(order=1))
lines!(ax, c, label="First-order UpwindBiased")

c = advect_tracer(WENO(order=11), closure=ScalarDiffusivity(κ=0.04))
lines!(ax, c, label="WENO(order=11) with κ=0.04")

Legend(fig[1, 2], ax, framevisible=false)
xlims!(ax, 1, 7)

save("advection_schemes.png", fig)
