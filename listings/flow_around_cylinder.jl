using Oceananigans
using Printf

config = :les
Re = 10^6
# config = :les

if config == :dns
    if Re == 1
        Ny = 256 
        Nx = 2Ny
    elseif Re == 10
        Ny = 256 
        Nx = 2Ny
    elseif Re == 100
        Ny = 512 
        Nx = 2Ny
    elseif Re == 1000
        Ny = 2^11
        Nx = 2Ny
    elseif Re == 10^4
        Ny = 2^12
        Nx = 2Ny
    elseif Re == 10^5
        Ny = 2^13
        Nx = 2Ny
    elseif Re == 10^6
        Ny = 3/2 * 2^13 |> Int
        Nx = 2Ny
    end

else
    Ny = 2^11
    Nx = 2Ny
    Re = Inf
end

prefix = "flow_around_cylinder_$(config)_Re$(Re)_Ny$(Ny)"

grid = RectilinearGrid(GPU(), size=(Nx, Ny), x=(-6, 42), y=(-12, 12), halo = (7, 7),
                       topology=(Periodic, Bounded, Flat))

cylinder(x, y) = (x^2 + y^2) ≤ 1
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(cylinder))

if config == :dns
    advection = Centered(order=2)
    closure = ScalarDiffusivity(ν=1/Re)

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

rate = 2
@inline mask(x, y, δ=6, x₀=42-δ) = max(zero(x), (x - x₀) / δ)
u_sponge = Relaxation(target=1; mask, rate)
v_sponge = Relaxation(target=0; mask, rate)
forcing = (u=u_sponge, v=v_sponge)

model = NonhydrostaticModel(; grid, closure,
                            advection, forcing, boundary_conditions, timestepper=:RungeKutta3)
uᵢ(x, y) = 1 + 1e-2 * randn()
set!(model, u=uᵢ)

Δx = minimum_xspacing(grid)
max_Δt = if config == :dns
    0.2 * Δx^2 * Re
else
    0.2 * Δx
end

simulation = Simulation(model, Δt=max_Δt, stop_time=200)
conjure_time_step_wizard!(simulation, cfl=0.5, IterationInterval(10); max_Δt)

progress(sim) = @info @sprintf("Iter: %d, time: %.2f", iteration(sim), time(sim))
add_callback!(simulation, progress, IterationInterval(100))

u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)
outputs = (; u, v, ζ)
simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(10),
                                                    filename = prefix * "_fields.jld2",
                                                    overwrite_existing = true,
                                                    with_halos = true)

# Drag computation
μ = XFaceField(grid)
set!(μ, mask)
r = u_sponge.rate
drag_force = Field(Integral(r * μ * (1 - u)))
compute!(drag_force)

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

