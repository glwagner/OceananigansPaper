using Oceananigans
using Oceananigans.Models: seawater_density
using SeawaterPolynomials: TEOS10EquationOfState

grid = RectilinearGrid(GPU(),
                       size = (2048, 2048),
                       x = (-0.25, 0.25),
                       z = (-0.5, 0.0),
                       topology = (Bounded, Flat, Bounded))

closure = ScalarDiffusivity(ν=1.15e-6, κ=(T=1e-7, S=1e-9))

equation_of_state = TEOS10EquationOfState(reference_density=1000)
buoyancy = SeawaterBuoyancy(; equation_of_state, gravitational_acceleration=9.81) 

model = NonhydrostaticModel(; grid, buoyancy, closure, tracers = (:T, :S), timestepper=:RungeKutta3)
                            
Tᵢ(x, z) = 1 + (7.55 - 1) * (z > -0.25)
Ξᵢ(x, z) = 1e-2 * randn()
set!(model, T=Tᵢ, S=0, u=Ξᵢ, v=Ξᵢ, w=Ξᵢ) 

Δx = minimum_xspacing(grid)
Δt = 0.2 * Δx^2 / closure.ν
simulation = Simulation(model; Δt, stop_time=40)
conjure_time_step_wizard!(simulation, cfl=0.5, max_Δt=Δt)

run!(simulation)

ρ = Field(seawater_density(model))
compute!(ρ)

using CairoMakie
heatmap!(axρ, ρ, colormap=:grays)
