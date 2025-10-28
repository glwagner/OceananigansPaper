using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials: TEOS10EquationOfState

H, L = 256meters, 1024meters
δ = L / 2
x, y, z = (-3L, 3L), (-L, L), (-H, 0)
Nz = 64

grid = RectilinearGrid(GPU(); size=(6Nz, 2Nz, Nz), halo=(6, 6, 6),
                       x, y, z, topology=(Bounded, Bounded, Bounded))

wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(wedge))

T₂ = 12.421hours
U₂ = 0.1 # m/s

@inline U(x, y, z, t, p) = p.U₂ * sin(2π * t / p.T₂)
@inline U(y, z, t, p) = U(zero(y), y, z, t, p)

open_bc = OpenBoundaryCondition(U; inflow_timescale = 2minutes,
                                   outflow_timescale = 2minutes,
                                   parameters=(; U₂, T₂),
                                   scheme = PerturbationAdvection())

u_bcs = FieldBoundaryConditions(east = open_bc, west = open_bc)

@inline ambient_temperature(x, z, t, H) = 12 + 4z/H
@inline ambient_temperature(x, y, z, t, H) = ambient_temperature(x, z, t, H)
ambient_temperature_bc = ValueBoundaryCondition(ambient_temperature; parameters = H)
T_bcs = FieldBoundaryConditions(east = ambient_temperature_bc, west = ambient_temperature_bc)

ambient_salinity_bc = ValueBoundaryCondition(32)
S_bcs = FieldBoundaryConditions(east = ambient_salinity_bc, west = ambient_salinity_bc)

model = NonhydrostaticModel(; grid, tracers = (:T, :S),
                              buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState()),
                              advection = WENO(order=9), coriolis = FPlane(latitude=47.5),
                              boundary_conditions = (; T=T_bcs, u = u_bcs, S = S_bcs))

Tᵢ(x, y, z) = ambient_temperature(x, y, z, 0, H)

set!(model, T=Tᵢ, S=32, u=U(0, 0, 0, 0, (; U₂, T₂)))

simulation = Simulation(model, Δt=5, stop_time=3days)
conjure_time_step_wizard!(simulation, cfl=0.7)

using Printf

wallclock = Ref(time_ns())
T₀ = time_ns()

function progress(sim)
    u, v, w = sim.model.velocities
    ΔT = 1e-9 * (time_ns() - wallclock[])
    ΔT₀ = 1e-9 * (time_ns() - T₀)
    msg = @sprintf("(%d) t: %s, Δt: %s, wall time since last print: %s, total wall time : %s",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(ΔT), prettytime(ΔT₀))
    @info msg
    wallclock[] = time_ns()
    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

prefix = "flow_past_headland_$Nz"
u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)

using Oceanostics: ErtelPotentialVorticity
using Oceananigans.BuoyancyFormulations: buoyancy
q = Field(ErtelPotentialVorticity(model, model.velocities..., buoyancy(model), model.coriolis))
outputs = merge(model.velocities, model.tracers, (; ζ, q))

xy_writer = JLD2Writer(model, outputs,
                       indices = (:, :, grid.Nz),
                       schedule = TimeInterval(10minutes),
                       filename = prefix * "_xy.jld2",
                       overwrite_existing = true)

xz_writer = JLD2Writer(model, outputs,
                       indices = (:, grid.Ny÷2, :),
                       schedule = TimeInterval(10minutes),
                       filename = prefix * "_xz.jld2",
                       overwrite_existing = true)

xyz_writer = JLD2Writer(model, outputs,
                        schedule = TimeInterval(20minutes),
                        filename = prefix * "_xyz.jld2",
                        overwrite_existing = true)

simulation.output_writers[:xy] = xy_writer
simulation.output_writers[:xz] = xz_writer
simulation.output_writers[:xyz] = xyz_writer

run!(simulation)
