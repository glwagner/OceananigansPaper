# Headland LES SYPD benchmark (paper config: 384×128×64, WENO(9), tidal forcing).
# Mirrors listings/flow_past_headland.jl exactly except that it uses the SYPD
# harness for short timed measurement instead of running the full 3 days.

using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials: TEOS10EquationOfState
using CUDA
using Printf

include(joinpath(@__DIR__, "sypd_harness.jl"))

H, L = 256meters, 1024meters
δ = L / 2
x, y, z = (-3L, 3L), (-L, L), (-H, 0)
Nz = parse(Int, get(ENV, "HEADLAND_NZ", "64"))

grid = RectilinearGrid(GPU(); size=(6Nz, 2Nz, Nz), halo=(6, 6, 6),
                       x, y, z, topology=(Bounded, Bounded, Bounded))

wedge(x, y) = -H * (1 + (y + abs(x)) / δ)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(wedge))

T₂ = 12.421hours
U₂ = 0.1

@inline U(x, y, z, t, p) = p.U₂ * sin(2π * t / p.T₂)
@inline U(y, z, t, p) = U(zero(y), y, z, t, p)

open_bc = OpenBoundaryCondition(U; parameters=(; U₂, T₂),
    scheme = PerturbationAdvection(inflow_timescale=2minutes, outflow_timescale=2minutes))

u_bcs = FieldBoundaryConditions(east=open_bc, west=open_bc)

@inline ambient_temperature(x, z, t, H) = 12 + 4z/H
@inline ambient_temperature(x, y, z, t, H) = ambient_temperature(x, z, t, H)
ambient_temperature_bc = ValueBoundaryCondition(ambient_temperature; parameters=H)
T_bcs = FieldBoundaryConditions(east=ambient_temperature_bc, west=ambient_temperature_bc)

ambient_salinity_bc = ValueBoundaryCondition(32)
S_bcs = FieldBoundaryConditions(east=ambient_salinity_bc, west=ambient_salinity_bc)

model = NonhydrostaticModel(grid; tracers=(:T, :S),
                            buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState()),
                            advection = WENO(order=9), coriolis=FPlane(latitude=47.5),
                            boundary_conditions = (; T=T_bcs, u=u_bcs, S=S_bcs))

Tᵢ(x, y, z) = ambient_temperature(x, y, z, 0, H)
set!(model, T=Tᵢ, S=32, u=U(0, 0, 0, 0, (; U₂, T₂)))

const Δt_seconds = parse(Float64, get(ENV, "HEADLAND_DT", "5"))
simulation = Simulation(model, Δt=Δt_seconds)
# Fixed Δt for the benchmark — adaptive wizard drives Δt to NaN once
# velocities develop. The listing's initial Δt=5 s is representative for
# the paper resolution + 0.1 m/s tidal forcing.

measure_sypd!(simulation; name="flow_past_headland",
              extra_config=(; Nz, weno_order=9, dt_seconds=Δt_seconds, U₂, T₂_seconds=Float64(T₂)))
