using Oceananigans
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization, VerticallyImplicitTimeDiscretization
using GLMakie

fig = Figure(size=(1000, 400))
axκ = Axis(fig[1, 1], xlabel="Tracer diffusivity", ylabel="z")
axc = Axis(fig[1, 2], xlabel="Tracer distribution at t=10", ylabel="z", yaxisposition=:right)

for (Δt, time_discretization) in [(10,   ExplicitTimeDiscretization()),
                                  (1e-4, ExplicitTimeDiscretization()),
                                  (10,   VerticallyImplicitTimeDiscretization())]

    grid = RectilinearGrid(size=200, z=(-5, 5), topology=(Flat, Flat, Bounded))

    κ(z, t) = exp(-z^2)
    closure = VerticalScalarDiffusivity(time_discretization; κ)

    model = HydrostaticFreeSurfaceModel(; grid, closure, tracers=:c,
                                        velocities=PrescribedVelocityFields())

    cᵢ(z) = z
    set!(model, c=cᵢ)

    simulation = Simulation(model; Δt, stop_time=10)
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

Legend(fig[1, 3], axc)
display(current_figure())

hidespines!(axκ, :t, :r)
hidespines!(axc, :t, :l)

save("vertically_implicit_diffusion.png", fig)

