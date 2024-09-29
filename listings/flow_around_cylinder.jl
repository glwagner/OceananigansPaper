using Oceananigans
using Oceananigans.Models.NonhydrostaticModels: ConjugateGradientPoissonSolver
using Printf
using CUDA

config = :les
u∞ = 1
r = 1/2
arch = GPU()
stop_time = 100

if config == :dns
    if Re <= 10
        Ny = 256 
        Nx = 2Ny
    elseif Re <= 100
        Ny = 512 
        Nx = 2Ny
    elseif Re <= 1000
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
    Ny = 512
    Nx = 2Ny
    Re = Inf
end

prefix = "flow_around_cylinder_$(config)_Re$(Re)_Ny$(Ny)"

x = (-3, 21) # 24
y = (-6, 6)  # 12

grid = RectilinearGrid(arch, size=(Nx, Ny); x, y, halo=(9, 9),
                       topology=(Periodic, Bounded, Flat))

cylinder(x, y) = (x^2 + y^2) ≤ r^2
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(cylinder))

@inline u_drag(x, y, t, u, v, Cᴰ) = - Cᴰ * u * sqrt(u^2 + v^2)
@inline v_drag(x, y, t, u, v, Cᴰ) = - Cᴰ * v * sqrt(u^2 + v^2)

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
    ℓ = 1e-4
    @show Cᴰ = (ϰ / log(x₁ / ℓ))^2
    u_drag_bc = FluxBoundaryCondition(u_drag, field_dependencies=(:u, :v), parameters=Cᴰ)
    v_drag_bc = FluxBoundaryCondition(v_drag, field_dependencies=(:u, :v), parameters=Cᴰ)
    u_bcs = FieldBoundaryConditions(immersed=u_drag_bc)
    v_bcs = FieldBoundaryConditions(immersed=v_drag_bc)
    boundary_conditions = (u=u_bcs, v=v_bcs)
end

rate = 1
x = xnodes(grid, Face())
@inline mask(x, y, δ=6, x₀=x[end]) = max(zero(x), 1 + (x - x₀) / δ)
u_sponge = Relaxation(target=1; mask, rate)
v_sponge = Relaxation(target=0; mask, rate)
forcing = (u=u_sponge, v=v_sponge)
#pressure_solver = ConjugateGradientPoissonSolver(grid, maxiter=100, reltol=eps(grid))
#pressure_solver = ConjugateGradientPoissonSolver(grid, maxiter=100)
pressure_solver = nothing

if isnothing(pressure_solver)
    prefix *= "_fft"
end

model = NonhydrostaticModel(; grid, pressure_solver, closure,
                            advection, forcing, boundary_conditions,
                            timestepper=:RungeKutta3)

@show model

uᵢ(x, y) = 1 + 1e-2 * randn()
set!(model, u=uᵢ)

Δx = minimum_xspacing(grid)
if config == :dns
    Δt = max_Δt = 0.2 * Δx^2 * Re
else
    Δt = 0.2 * Δx
    max_Δt = Inf
end

simulation = Simulation(model; Δt, stop_time)
conjure_time_step_wizard!(simulation, cfl=0.7, IterationInterval(3); max_Δt)

u, v, w = model.velocities
d = ∂x(u) + ∂y(v)

# Drag computation
μ = XFaceField(grid)
set!(μ, mask)
ω = u_sponge.rate
drag_force = Field(Integral(ω * μ * (u∞ - u)))
compute!(drag_force)

function progress(sim)
    if pressure_solver isa ConjugateGradientPoissonSolver
        pressure_iters = iteration(pressure_solver)
    else
        pressure_iters = 0
    end

    compute!(drag_force)
    D = CUDA.@allowscalar first(drag_force)
    cᴰ = D / (u∞ * r) 
    vmax = maximum(model.velocities.v)
    dmax = maximum(abs, d)

    @info @sprintf("Iter: %d, time: %.2f, Δt: %.4f, Poisson iters: %d, max d: %.2e, max v: %.2e, cᴰ = %0.2f",
                   iteration(sim), time(sim), sim.Δt, pressure_iters, dmax, vmax, cᴰ)

    return nothing
end

add_callback!(simulation, progress, IterationInterval(1))

ζ = ∂x(v) - ∂y(u)
outputs = (; u, v, ζ, d)

simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(10),
                                                    filename = prefix * "_fields.jld2",
                                                    overwrite_existing = true,
                                                    with_halos = true)

simulation.output_writers[:drag] = JLD2OutputWriter(model, (; drag_force),
                                                    schedule = TimeInterval(0.1),
                                                    filename = prefix * "_drag.jld2",
                                                    overwrite_existing = true,
                                                    with_halos = true)
 
run!(simulation)

#=
using GLMakie
ζ = Field(ζ)
compute!(ζ)

fig = Figure(size=(1600, 500))
ax = Axis(fig[1, 1], aspect=3)
heatmap!(ax, ζ, colormap=:balance, nan_color=:gray)
display(fig)
=#

