using Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: PartialCellBottom, mask_immersed_field!
using Oceananigans.Models.NonhydrostaticModels: ImmersedPoissonSolver, DiagonallyDominantPreconditioner
using Printf

Lx = 512 # domain depth (m)
Lz = 512 # domain depth (m)
h₀ = 10  # characteristic mountain height
δ  = Lx / 100

grid = RectilinearGrid(size = (256, 256), halo = (4, 4),
                       x = (-Lx/2, Lx/2), z = (-Lz, 0),
                       topology = (Periodic, Flat, Bounded))

using Random: seed!
seed!(123)
seamounts = 42
W = Lx - 2δ
positions = W  .* (rand(seamounts) .- 1/2) # ∈ [-Lx, Lx]
heights   = h₀ .* (rand(seamounts) .+ 1) # ∈ [h₀, 2h₀]

function bottom(x)
    z = - Lz
    for (h, ξ) in zip(heights, positions)
        z += h * exp(-(x - ξ)^2 / 2δ^2)
    end
    return z
end

grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))

@inline tidal_forcing(x, z, t, p) = p.U₂ * 2π/p.T₂ * sin(2π / p.T₂ * t)
u_forcing = Forcing(tidal_forcing, parameters=(; U₂=0.1, T₂=12.421hours))

b_immersed_bc = FluxBoundaryCondition(5e-8)
b_bcs = FieldBoundaryConditions(immersed=b_immersed_bc)
boundary_conditions = (; b=b_bcs)

pressure_solver = ImmersedPoissonSolver(grid)#, maxiter=3)
#pressure_solver = ImmersedPoissonSolver(grid, preconditioner=DiagonallyDominantPreconditioner())
#pressure_solver = ImmersedPoissonSolver(grid, preconditioner=nothing)

@show pressure_solver

model = NonhydrostaticModel(; grid, boundary_conditions,
                            pressure_solver,
                            tracers = :b,
                            timestepper = :RungeKutta3,
                            buoyancy = BuoyancyTracer(),
                            advection = WENO())
                            #forcing = (; u = u_forcing))

bᵢ(x, z) = 1e-6 * z
set!(model, b=bᵢ)

T₂ = u_forcing.parameters.T₂

simulation = Simulation(model; Δt=1.0, stop_time=12hours)
conjure_time_step_wizard!(simulation, cfl=0.7)

u, v, w = model.velocities
δ = ∂x(u) + ∂z(w)

function progress(sim)
    cg_iter = model.pressure_solver.pcg_solver.iteration
    δmax = maximum(δ)
    msg = @sprintf("Sim iter: %d, cg iter: %d, max δ: %.2e, time: %s, Δt: %s",
                   iteration(sim), cg_iter, δmax, prettytime(sim), prettytime(sim.Δt))
    @info msg
    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)

using GLMakie

δf = Field(δ)
compute!(δf)

fig = Figure(size=(1400, 500))

axw = Axis(fig[1, 1], xlabel="x (m)", ylabel="z (m)", aspect=1)
axb = Axis(fig[1, 2], xlabel="x (m)", ylabel="z (m)", aspect=1)
axd = Axis(fig[1, 3], xlabel="x (m)", ylabel="z (m)", aspect=1)

wlim = 4e-2 #maximum(abs, w) / 2
b = model.tracers.b

heatmap!(axw, w, colormap=:balance, nan_color=:gray, colorrange=(-wlim, wlim))
heatmap!(axb, b, colormap=:magma, nan_color=:gray)
heatmap!(axd, δf, colormap=:balance, nan_color=:gray)

display(fig)

