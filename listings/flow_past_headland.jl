using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials

H, L, δ = 256meters, 1024meters, 512meters
x, y, z = (-3L, 3L), (-L, L), (-H, 0)
Nz = 64

grid = RectilinearGrid(GPU(); size=(6Nz, 2Nz, Nz), halo=(6, 6, 6),
                       x, y, z, topology=(Periodic, Bounded, Bounded))

wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(wedge))

@show grid

T₂ = 12.421hours
U₂ = 0.1 # m/s
@inline Fu(x, y, z, t, p) = 2π * p.U₂ / p.T₂ * sin(2π * t / p.T₂)
forcing = (; u=Forcing(Fu, parameters=(; U₂, T₂)))

equation_of_state = SeawaterPolynomials.TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)

Q = 250.0  # W m⁻², surface _heat_ flux
ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater
Jᵀ = Q / (ρₒ * cᴾ) # K m s⁻¹, surface _temperature_ flux
surface_temperature_flux(x, y, t, p) = p.Jᵀ * (1 - exp(-t^4 / (p.time_delay^4)))
temperature_top_bc = FluxBoundaryCondition(surface_temperature_flux, parameters = (; Jᵀ, time_delay = 1day))
temperature_bcs = FieldBoundaryConditions(top = temperature_top_bc)

model = NonhydrostaticModel(; grid, tracers = (:T, :S), buoyancy, forcing,
                            advection = WENO(order=9), coriolis = FPlane(latitude=47.5),
                            boundary_conditions = (; T=temperature_bcs))

Tᵢ(x, y, z) = 12 + 4z / H
set!(model, T=Tᵢ, S=32, u=U₂/2)

simulation = Simulation(model, Δt=5, stop_time=3days)
conjure_time_step_wizard!(simulation, cfl=0.7)

using Printf

wallclock = Ref(time_ns())
T₀ = time_ns()

function progress(sim)
    u, v, w = sim.model.velocities
    ΔT = 1e-9 * (time_ns() - wallclock[])
    ΔT₀ = 1e-9 * (time_ns() - T₀)
    msg = @sprintf("(%d) t: %s, Δt: %s, wall time / time step: %s, total wall time : %s",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(ΔT), prettytime(ΔT₀))
    @info msg
    wallclock[] = time_ns()
    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

prefix = "flow_past_headland_$Nz"
u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)
s = @at (Center, Center, Center) sqrt(u^2 + v^2)

using Oceanostics: ErtelPotentialVorticity
q = ErtelPotentialVorticity(model, model.velocities..., model.tracers.T, model.coriolis)
outputs = merge(model.velocities, model.tracers, (; ζ, s, q))

xy_writer = JLD2OutputWriter(model, outputs,
                             indices = (:, :, grid.Nz),
                             schedule = TimeInterval(10minutes),
                             filename = prefix * "_xy.jld2",
                             overwrite_existing = true)

xz_writer = JLD2OutputWriter(model, outputs,
                             indices = (:, grid.Ny÷2, :),
                             schedule = TimeInterval(10minutes),
                             filename = prefix * "_xz.jld2",
                             overwrite_existing = true)
xyz_writer = JLD2OutputWriter(model, outputs,
                              schedule = TimeInterval(12hours),
                              filename = prefix * "_xyz.jld2",
                              overwrite_existing = true)

simulation.output_writers[:xy] = xy_writer
simulation.output_writers[:xz] = xz_writer
simulation.output_writers[:xyz] = xyz_writer

run!(simulation)
