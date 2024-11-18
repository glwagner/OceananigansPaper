using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures: CATKEVerticalDiffusivity
using ClimaOcean
using SeawaterPolynomials: TEOS10EquationOfState
using Printf

arch = CPU()
Nx = Ny = 128
Nz = 32

# Geometry parameters
H∞ = 3kilometers
Wˢ = 150kilometers
Hˢ = 500meters
y₀ = -3Wˢ/2 + 50kilometers

x = (0, 3Wˢ)
y = (-3Wˢ/2, 3Wˢ/2)
z = exponential_z_faces(; Nz, depth=H∞)

grid = RectilinearGrid(arch; size=(Nx, Ny, Nz), halo=(7, 7, 7),
                       x, y, z, topology=(Periodic, Bounded, Bounded))
                       
bottom_height(x, y) = -Hˢ - (H∞ - Hˢ) * (1 + tanh(y / Wˢ)) / 2
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

@inline τx(x, y, t, p) = - p.τ₀ * sin(π * (y - p.y₀) / (p.yᴱ - p.y₀)) * (y > p.y₀)
u_top_bc = FluxBoundaryCondition(τx, parameters=(τ₀=1e-4, y₀=y₀, yᴱ=y[2]))

@inline Jˢ(x, y, t, p) = p.Jˢ * (y < p.y₀)
S_top_bc = FluxBoundaryCondition(Jˢ, parameters=(Jˢ=-3e-3, y₀=y₀))

u_bcs = FieldBoundaryConditions(top=u_top_bc)
S_bcs = FieldBoundaryConditions(top=S_top_bc)

equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

model = HydrostaticFreeSurfaceModel(; grid, buoyancy,
                                    tracers = (:T, :S, :e),
                                    closure = CATKEVerticalDiffusivity(),
                                    momentum_advection = WENOVectorInvariant(),
                                    tracer_advection = WENO(order=9),
                                    coriolis = BetaPlane(latitude=-60),
                                    boundary_conditions = (u=u_bcs, S=S_bcs))

w, h = 50kilometers, 1000meters
Tᵢ(x, y, z) = 5 * exp(z / h) + tanh(y / w) + 1e-2 * randn()
Sᵢ(x, y, z) = 32 + 0.5 * exp(z / h) + 1e-3 * randn()
set!(model, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=1minute, stop_iteration=100)

wallclock = Ref(time_ns())

function progress(sim)
    T = sim.model.tracers.T
    S = sim.model.tracers.S
    elapsed = 1e-9 * (time_ns() - wallclock[])
    msg = @sprintf("%d | t: %s, wall time: %s, extrema(T): (%.2f, %.2f), extrema(S): (%.2f, %.2f)",
                   iteration(sim), prettytime(sim), prettytime(elapsed),
                   minimum(T), maximum(T), minimum(S), maximum(S))

    u, v, w = sim.model.velocities
    msg *= @sprintf("max|u|: (%.2e, %.2e, %.2e)",
                    maximum(abs, u), maximum(abs, v), maximum(abs, w))
                    
    @info msg
                   
    wallclock[] = time_ns()
    return nothing
end

xy = JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                      filename = "antarctic_slope_current.jld2",
                      schedule = TimeInterval(6hours),
                      indices = (:, :, Nz),
                      overwrite_existing = true)

simulation.output_writers[:xy] = xy

run!(simulation)

