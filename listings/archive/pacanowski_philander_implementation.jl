using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization

grid = RectilinearGrid(size=50, z=(-200, 0), topology=(Flat, Flat, Bounded))

ν = Field{Center, Center, Face}(grid)
κ = Field{Center, Center, Face}(grid)
vitd = VerticallyImplicitTimeDiscretization()
closure = VerticalScalarDiffusivity(vitd; ν, κ)

b_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(1e-7))
u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(-5e-4))

model = HydrostaticFreeSurfaceModel(; grid, closure, tracers=:b,
                                    boundary_conditions = (; b=b_bcs, u=u_bcs),
                                    buoyancy = BuoyancyTracer())

N² = 1e-5
bᵢ(z) = N² * z
set!(model, b=bᵢ)

simulation = Simulation(model, Δt=10minute, stop_time=2day)

b = model.tracers.b
u, v, w = model.velocities
Ri = Field(∂z(b) / (∂z(u)^2 + ∂z(v)^2))

function compute_diffusivities!(sim)
    ν₀ = 1e-2
     n = 2
     α = 5

    compute!(Ri)

    Riᵖ = parent(Ri)
    Riᵖ[isnan.(Riᵖ)] .= 0

    @. ν = ν₀ / (1 + α * Ri)^n
    @. κ = ν₀ / (1 + α * Ri)^(n+1)

    return nothing
end

add_callback!(simulation, compute_diffusivities!) 

run!(simulation)

using GLMakie

set_theme!(Theme(linewidth=3, linealpha=0.6))

fig = Figure(size=(800, 400))
axb = Axis(fig[1, 1], xlabel="Buoyancy", ylabel="z (m)")
axu = Axis(fig[1, 2], xlabel="x-velocity")
axκ = Axis(fig[1, 3], xlabel="Diffusivities", ylabel="z (m)", yaxisposition=:right)

lines!(axb, b)
lines!(axu, u)
lines!(axκ, ν, label="ν")
lines!(axκ, κ, label="κ")
axislegend(axκ, position=:lt)

hidespines!(axb, :t, :r)
hidespines!(axu, :t, :r, :l)
hidespines!(axκ, :t, :l)

hideydecorations!(axu, grid=true)

display(fig)

save("pacanowski_philander.png", fig)
