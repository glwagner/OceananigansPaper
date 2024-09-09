using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures:
    CATKEVerticalDiffusivity,
    TKEDissipationVerticalDiffusivity,
    RiBasedVerticalDiffusivity

function run_vertical_mixing_simulation(closure)
    Jb = 1e-7
    τx = -5e-4
    N² = 1e-5

    grid = RectilinearGrid(size=50, z=(-200, 0), topology=(Flat, Flat, Bounded))

    b_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jb))
    u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(τx))

    if closure isa CATKEVerticalDiffusivity
        tracers = (:b, :e)
    elseif closure isa TKEDissipationVerticalDiffusivity
        tracers = (:b, :e, :ϵ)
    else
        tracers = :b
    end

    model = HydrostaticFreeSurfaceModel(; grid, closure, tracers,
                                        boundary_conditions = (; b=b_bcs, u=u_bcs),
                                        buoyancy = BuoyancyTracer())

    bᵢ(z) = N² * z
    set!(model, b=bᵢ)

    simulation = Simulation(model, Δt=1minute, stop_time=24hours)
    run!(simulation)

    return simulation
end

closures = Dict("CATKE" => CATKEVerticalDiffusivity(),
                "k-epsilon" => TKEDissipationVerticalDiffusivity(),
                "Ri-based" => RiBasedVerticalDiffusivity())

using GLMakie

fig = Figure(size=(800, 600))
axb = Axis(fig[1, 1], xlabel="Buoyancy (m s⁻²)", ylabel="z (m)")
axu = Axis(fig[1, 2], xlabel="x-velocity (m s⁻¹)", ylabel="z (m)")
axκ = Axis(fig[1, 3], xlabel="Tracer diffusivity (m² s⁻¹)", ylabel="z (m)")

for (name, closure) in closures
    sim = run_vertical_mixing_simulation(closure)
    lines!(axb, sim.model.tracers.b, label=name)
    lines!(axu, sim.model.velocities.u, label=name)
    lines!(axκ, sim.model.diffusivity_fields.κc, label=name)
end

axislegend(axb, position=:lt)
display(fig)

