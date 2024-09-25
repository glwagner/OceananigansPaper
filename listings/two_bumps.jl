using Oceananigans
using Oceananigans.Units
using Oceananigans.Solvers: solve!
using Oceananigans.ImmersedBoundaries: PartialCellBottom, mask_immersed_field!
using Oceananigans.Models.NonhydrostaticModels: PressureSolver, ImmersedPoissonSolver, DiagonallyDominantPreconditioner
using Printf

H  = 2kilometers  # domain depth (m)
h₁ = 500          # first mountain height (m)
h₂ = 1000         # second mountain height (m)
δ  = 1kilometers # mountain width (m)
L  = 10kilometers

underlying_grid = RectilinearGrid(size = (128, 64),
                                  x = (-L, L),
                                  z = (-H, 0),
                                  halo = (4, 4),
                                  topology = (Periodic, Flat, Bounded))

N² = 1e-5
f = 1e-4
@show λ = π * sqrt(N²) * H / f # first mode wavelength
x₁ = + L/4 #2λ
x₂ = - L/4 #2λ

bottom(x) = -H + h₁ * exp(-(x - x₁)^2 / 2δ^2) +
                 h₂ * exp(-(x - x₂)^2 / 2δ^2)

#grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom))
grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom))

@inline tidal_forcing(x, z, t, p) = p.U₂ * 2π/p.T₂ * sin(2π / p.T₂ * t)
u_forcing = Forcing(tidal_forcing, parameters=(; U₂=0.1, T₂=12.421hours))

pressure_solver = ImmersedPoissonSolver(grid) #, maxiter=10) #, maxiter=Inf)
#pressure_solver = ImmersedPoissonSolver(grid, preconditioner=DiagonallyDominantPreconditioner())
#pressure_solver = ImmersedPoissonSolver(grid, preconditioner=nothing)

model = NonhydrostaticModel(; grid,
                            pressure_solver,
                            timestepper = :RungeKutta3,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            coriolis = FPlane(; f),
                            advection = WENO(),
                            forcing = (; u=u_forcing))

@show model.pressure_solver

bᵢ(x, z) = 1e-5 * z
set!(model, b=bᵢ)
T₂ = 12.421hours
simulation = Simulation(model; Δt=5minutes, stop_time=2T₂)

u, v, w = model.velocities
δ = ∂x(u) + ∂z(w)

function progress(sim)
    i = model.pressure_solver.pcg_solver.iteration
    δmax = maximum(δ)
    msg = @sprintf("Sim iter: %d, IPS iter: %d, max δ: %.2e", iteration(sim), i, δmax)
    @info msg
    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

run!(simulation)

using GLMakie

δf = Field(δ)
compute!(δf)

fig = Figure(size=(1200, 600))

axw = Axis(fig[1, 1], xlabel="x (m)", ylabel="z (m)", aspect=6)
axd = Axis(fig[2, 1], xlabel="x (m)", ylabel="z (m)", aspect=6)

heatmap!(axw, w, colormap=:balance, nan_color=:gray)
heatmap!(axd, δf, colormap=:balance, nan_color=:gray)

display(fig)

