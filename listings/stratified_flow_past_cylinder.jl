using Oceananigans
using Printf

# config = :les
config = :dns
Re = 100
Fr = 0.1
const N² = 1 / Fr^2
@inline bᵢ(x, z, t=0) = N² * z

if config == :dns
    if Re == 1
        Nz = 256 
        Nx = 3Nz
    elseif Re == 10
        Nz = 256 
        Nx = 3Nz
    elseif Re == 100
        Nz = 512 
        Nx = 3Nz
    elseif Re == 1000
        Nz = 2^11
        Nx = 3Nz
    end
elseif config == :les
    Nz = 512
    Nx = 3Nz
    Re = Inf
end

prefix = "stratified_flow_over_cylinder_$(config)_Re$(Re)_Nz$(Nz)"

grid = RectilinearGrid(GPU(),
                       size=(Nx, Nz), x=(-3, 12), z=(0, 5), halo = (7, 7),
                       topology=(Periodic, Flat, Bounded))

cylinder(x, z) = (x^2 + z^2) ≤ 1/2
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(cylinder))

if config == :dns
    advection = Centered(order=2)
    closure = ScalarDiffusivity(ν=1/Re, κ=1/Re)
    no_slip = ValueBoundaryCondition(0)
    velocity_bcs = FieldBoundaryConditions(immersed=no_slip)
    boundary_conditions = (u=velocity_bcs, v=velocity_bcs)

elseif config == :les
    advection = WENO(order=9)
    closure = nothing

    x₁ = minimum_xspacing(grid) / 2
    ϰ = 0.4
    ℓ = 1e-1
    Cᴰ = (ϰ / log(x₁ / ℓ))^2
    u_drag_bc = FluxBoundaryCondition(u_drag, field_dependencies=(:u, :v), parameters=Cᴰ)
    v_drag_bc = FluxBoundaryCondition(v_drag, field_dependencies=(:u, :v), parameters=Cᴰ)
    u_bcs = FieldBoundaryConditions(immersed=u_drag_bc)
    v_bcs = FieldBoundaryConditions(immersed=v_drag_bc)
    boundary_conditions = (u=u_bcs, v=v_bcs)
end


@inline mask(x, z, δ=2, x₀=10) = max(zero(x), (x - x₀) / δ)
u_sponge = Relaxation(rate=10; target=1, mask)
w_sponge = Relaxation(rate=10; target=0, mask)
b_sponge = Relaxation(rate=10; target=bᵢ, mask)
forcing = (u=u_sponge, w=w_sponge, b=b_sponge)

model = NonhydrostaticModel(; grid, closure, advection, forcing,
                            tracers = :b, buoyancy = BuoyancyTracer(),
                            boundary_conditions, timestepper=:RungeKutta3)

uᵢ(x, z) = 1 + 1e-2 * randn()
set!(model, u=uᵢ, b=bᵢ)

Δx = minimum_xspacing(grid)
max_Δt = if config == :dns
    0.2 * Δx^2 * Re
else
    0.2 * Δx
end

simulation = Simulation(model, Δt=max_Δt, stop_time=100)
conjure_time_step_wizard!(simulation, cfl=0.5, IterationInterval(10); max_Δt)

progress(sim) = @info @sprintf("Iter: %d, time: %.2f, max|w|: %.2e", iteration(sim), time(sim),
                               maximum(abs, interior(sim.model.velocities.w)))

add_callback!(simulation, progress, IterationInterval(100))

u, v, w = model.velocities
b = model.tracers.b
ξ = ∂z(u) - ∂x(w)
outputs = (; u, w, b, ξ)
simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(10),
                                                    filename = prefix * "_fields.jld2",
                                                    overwrite_existing = true,
                                                    with_halos = true)

# Drag computation
μ = XFaceField(grid)
set!(μ, mask)
r = u_sponge.rate # relaxation rate
u★ = 1
drag = Integral(r * μ * (u★ - u))

simulation.output_writers[:drag] = JLD2OutputWriter(model, (; drag),
                                                    schedule = TimeInterval(0.1),
                                                    filename = prefix * "_drag.jld2",
                                                    overwrite_existing = true,
                                                    with_halos = true)
 
run!(simulation)

#=
using GLMakie
u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u))
compute!(ζ)

heatmap(ζ, nan_color=:gray)
display(current_figure())
=#

