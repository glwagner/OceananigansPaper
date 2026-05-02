# Eady LES SYPD benchmark (paper config: 1000×1000×64, WENO(9)).
# Run on a single H100. Builds the model at paper resolution and measures SYPD
# over a short timed window via benchmarks/sypd_harness.jl. Writes JSON to
# benchmarks/runs/eady_les_<jobid>.json.

using Oceananigans
using Oceananigans.Units
using CUDA
using Printf

include(joinpath(@__DIR__, "sypd_harness.jl"))

const Nx = parse(Int, get(ENV, "EADY_NX", "1000"))
const Ny = parse(Int, get(ENV, "EADY_NY", "1000"))
const Nz = parse(Int, get(ENV, "EADY_NZ", "64"))
const Lx = parse(Float64, get(ENV, "EADY_LX", "4000"))
const Lz = parse(Float64, get(ENV, "EADY_LZ", "128"))
const weno_order = parse(Int, get(ENV, "EADY_WENO", "9"))

@info @sprintf("Eady benchmark config: %d×%d×%d, %g×%g×%g m, WENO(%d)", Nx, Ny, Nz, Lx, Lx, Lz, weno_order)

arch = GPU()
grid = RectilinearGrid(arch,
                       size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (0, Lx),
                       z = (-Lz, 0),
                       topology = (Periodic, Periodic, Bounded))

f = 1e-4
N² = 1e-6
Ri = 1
Λ = sqrt(N² / Ri)
parameters = (; Λ, f)

@inline uᵇᵍ(x, y, z, t, p) = + p.Λ * z
@inline bᵇᵍ(x, y, z, t, p) = - p.f * p.Λ * y
background_fields = (u = BackgroundField(uᵇᵍ; parameters),
                     b = BackgroundField(bᵇᵍ; parameters))

model = NonhydrostaticModel(grid; background_fields,
                            coriolis = FPlane(; f),
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            advection = WENO(order=weno_order))

Δz = minimum_zspacing(grid)
bᵢ(x, y, z) = N² * z + 1e-2 * N² * Δz * (2rand() - 1)
set!(model, b=bᵢ)

const Δt_seconds = parse(Float64, get(ENV, "EADY_DT", "30"))
simulation = Simulation(model, Δt=Δt_seconds)
# Fixed Δt for the benchmark — the wizard's default max_Δt=Inf lets the
# initial near-quiescent flow blow up Δt to NaN. A representative Δt for
# the paper's 4 m/2 m resolution + ~0.1 m/s velocities is ~30 s.

measure_sypd!(simulation; name="eady_les",
              extra_config=(; Nx, Ny, Nz, Lx, Lz, weno_order, dt_seconds=Δt_seconds))
