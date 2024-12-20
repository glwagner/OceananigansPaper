using Oceananigans
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, VerticallyImplicitTimeDiscretization
using GLMakie

fig = Figure(size=(1000, 400))
axκ = Axis(fig[1, 1], xlabel="Tracer diffusivity", ylabel="z")
axc = Axis(fig[1, 2], xlabel="Tracer distribution", ylabel="z", yaxisposition=:right)

for (Δt, time_discretization) in [(0.5,  ExplicitTimeDiscretization()),
                                  (1e-4, ExplicitTimeDiscretization()),
                                  (0.5,  VerticallyImplicitTimeDiscretization())]

    grid = RectilinearGrid(size=20, z=(-2, 2), topology=(Flat, Flat, Bounded))

    κ(z, t) = abs(z) < 1
    closure = VerticalScalarDiffusivity(time_discretization; κ)

    model = HydrostaticFreeSurfaceModel(; grid, closure, tracers=:c,
                                        velocities=PrescribedVelocityFields())

    cᵢ(z) = z
    set!(model, c=cᵢ)

    simulation = Simulation(model; Δt, stop_time=0.5)
    run!(simulation)
    label = string(summary(time_discretization), ", Δt=", Δt)
    lines!(axc, model.tracers.c; label)

    if time_discretization isa VerticallyImplicitTimeDiscretization # do this once
        z = znodes(model.tracers.c)
        lines!(axκ, κ.(z, 0), z, color=:black)

        z = znodes(model.tracers.c)
        lines!(axc, cᵢ.(z), z, linestyle=:dash, color=(:black, 0.5), label="Tracer distribution at t=0")
    end
end

Legend(fig[1, 3], axc, framevisible=false)
display(current_figure())

hidespines!(axκ, :t, :r)
hidespines!(axc, :t, :l)

xlims!(axc, -2, 2)
ylims!(axc, -2, 2)
ylims!(axκ, -2, 2)

save("vertically_implicit_diffusion.png", fig)

