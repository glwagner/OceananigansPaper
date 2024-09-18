using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures:
    CATKEVerticalDiffusivity,
    TKEDissipationVerticalDiffusivity,
    RiBasedVerticalDiffusivity

function run_vertical_mixing_simulation(closure)
    τx = -5e-4
    N² = 1e-5

    grid = RectilinearGrid(size=50, z=(-200, 0), topology=(Flat, Flat, Bounded))

    u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(τx))

    if closure isa CATKEVerticalDiffusivity
        tracers = (:b, :e)
    elseif closure isa TKEDissipationVerticalDiffusivity
        tracers = (:b, :e, :ϵ)
    else
        tracers = :b
    end

    model = HydrostaticFreeSurfaceModel(; grid, closure, tracers,
                                        boundary_conditions = (; u=u_bcs),
                                        buoyancy = BuoyancyTracer())

    bᵢ(z) = N² * z
    set!(model, b=bᵢ)

    simulation = Simulation(model, Δt=10, stop_time=24hours)
    run!(simulation)

    return simulation
end

closures = Dict("CATKE" => CATKEVerticalDiffusivity(),
                "k-epsilon" => TKEDissipationVerticalDiffusivity(),
                "Ri-based" => RiBasedVerticalDiffusivity())

using GLMakie

set_theme!(Theme(linewidth=3, linealpha=0.6))

fig = Figure(size=(800, 400))
axb = Axis(fig[1, 1], xlabel="Buoyancy (m s⁻²)", ylabel="z (m)")
axu = Axis(fig[1, 2], xlabel="x-velocity (m s⁻¹)", ylabel="z (m)")
axκ = Axis(fig[1, 3], xlabel="Tracer diffusivity (m² s⁻¹)", ylabel="z (m)", yaxisposition=:right)

for (name, closure) in closures
    sim = run_vertical_mixing_simulation(closure)
    lines!(axb, sim.model.tracers.b, label=name)
    lines!(axu, sim.model.velocities.u, label=name)
    lines!(axκ, sim.model.diffusivity_fields.κc, label=name)
end

axislegend(axb, position=:lt)

hidespines!(axb, :t, :r)
hidespines!(axu, :t, :r, :l)
hidespines!(axκ, :t, :l)

hideydecorations!(axu, grid=false)

display(fig)
save("vertical_mixing_parameterizations.png", fig)

