# Eady LES SYPD benchmark (paper production config: 1024×1024×64, WENO(9)).
# Adaptive time-stepping via `conjure_time_step_wizard!`. The Eady spinup starts
# nearly quiescent (random b perturbation only), so |𝐮| → 0 at t=0 and any
# pre-equilibrated SYPD reading would be meaningless: Δt is wherever the wizard
# happens to be on its way up, not where the equilibrated front-MLI cascade
# pins it. We therefore integrate forward by `EADY_SPINUP_DAYS` (default 5)
# of physical time before handing off to `measure_sypd!`, so the mean Δt the
# harness records is the equilibrated one.
#
# Single H100; writes benchmarks/runs/eady_les_<jobid>.json.

using Oceananigans
using Oceananigans.Units
using CUDA
using Printf

include(joinpath(@__DIR__, "sypd_harness.jl"))

const Nx = parse(Int,     get(ENV, "EADY_NX",          "1024"))
const Ny = parse(Int,     get(ENV, "EADY_NY",          "1024"))
const Nz = parse(Int,     get(ENV, "EADY_NZ",          "64"))
const Lx = parse(Float64, get(ENV, "EADY_LX",          "4096"))
const Lz = parse(Float64, get(ENV, "EADY_LZ",          "128"))
const weno_order  = parse(Int,     get(ENV, "EADY_WENO",        "9"))
const initial_dt  = parse(Float64, get(ENV, "EADY_DT0",         "30"))   # seed Δt (s)
const max_dt      = parse(Float64, get(ENV, "EADY_MAX_DT",      "120"))  # wizard ceiling (s); finite to keep the near-quiescent startup from blowing Δt → NaN
const cfl         = parse(Float64, get(ENV, "EADY_CFL",         "0.7"))
const spinup_days = parse(Float64, get(ENV, "EADY_SPINUP_DAYS", "5"))

@info @sprintf("Eady benchmark config: %d×%d×%d on %g×%g×%g m, WENO(%d), spinup %g days, CFL=%g, Δt₀=%g s, Δt_max=%g s",
               Nx, Ny, Nz, Lx, Lx, Lz, weno_order, spinup_days, cfl, initial_dt, max_dt)

arch = GPU()
grid = RectilinearGrid(arch,
                       size = (Nx, Ny, Nz),
                       x = (0, Lx),
                       y = (0, Lx),
                       z = (-Lz, 0),
                       topology = (Periodic, Periodic, Bounded))

f  = 1e-4
N² = 1e-6
Ri = 1
Λ  = sqrt(N² / Ri)
parameters = (; Λ, f)

@inline uᵇᵍ(x, y, z, t, p) = + p.Λ * z
@inline bᵇᵍ(x, y, z, t, p) = - p.f * p.Λ * y
background_fields = (u = BackgroundField(uᵇᵍ; parameters),
                     b = BackgroundField(bᵇᵍ; parameters))

model = NonhydrostaticModel(grid; background_fields,
                            coriolis  = FPlane(; f),
                            tracers   = :b,
                            buoyancy  = BuoyancyTracer(),
                            advection = WENO(order=weno_order))

Δz = minimum_zspacing(grid)
bᵢ(x, y, z) = N² * z + 1e-2 * N² * Δz * (2rand() - 1)
set!(model, b=bᵢ)

simulation = Simulation(model, Δt=initial_dt)
conjure_time_step_wizard!(simulation; cfl, max_Δt=max_dt)

# Progress callback. Reports an *instantaneous* SYPD over each 100-iter window
# so we can watch the wizard's Δt — and therefore SYPD — equilibrate during
# spinup rather than waiting until measure_sypd!'s timed window at the very
# end to find out whether it stabilized.
wallclock      = Ref(time_ns())
prev_sim_time  = Ref(simulation.model.clock.time)
prev_iteration = Ref(simulation.model.clock.iteration)
function progress(sim)
    now_wall   = time_ns()
    Δwall      = 1e-9 * (now_wall - wallclock[])
    Δsim_sec   = sim.model.clock.time      - prev_sim_time[]
    Δiter      = sim.model.clock.iteration - prev_iteration[]
    sypd_inst  = Δwall > 0 ? Δsim_sec / (Δwall * 365.25) : NaN
    sec_per_it = Δiter > 0 ? Δwall / Δiter : NaN
    u, v, w = sim.model.velocities
    @info @sprintf("iter %d, t = %s, Δt = %s, wall = %.2f s, sec/iter = %.3f, SYPD ≈ %.4f, max|u,v,w| = (%.3f, %.3f, %.3f) m/s",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt),
                   Δwall, sec_per_it, sypd_inst,
                   maximum(abs, u), maximum(abs, v), maximum(abs, w))
    wallclock[]      = now_wall
    prev_sim_time[]  = sim.model.clock.time
    prev_iteration[] = sim.model.clock.iteration
    return nothing
end
add_callback!(simulation, progress, IterationInterval(100))

# === Spinup phase: integrate by physical time until the wizard equilibrates ===
@info @sprintf("[spinup] integrating to t = %g days with adaptive Δt …", spinup_days)
empty!(simulation.output_writers)
empty!(simulation.diagnostics)
simulation.stop_time = spinup_days * 86400
simulation.stop_iteration = Inf
spinup_wall_start = time()
run!(simulation)
spinup_wall = time() - spinup_wall_start
spinup_iters = simulation.model.clock.iteration
spinup_final_dt = float(simulation.Δt)

u, v, w = simulation.model.velocities
isfinite(maximum(abs, u)) && isfinite(maximum(abs, v)) && isfinite(maximum(abs, w)) ||
    error("[spinup] NaN/Inf in velocities after $(spinup_days)-day spinup — pick a smaller EADY_MAX_DT.")

@info @sprintf("[spinup] done: %d iters in %.1f s wall, Δt now %s, max|u,v,w| = (%.3f, %.3f, %.3f) m/s",
               spinup_iters, spinup_wall, prettytime(spinup_final_dt),
               maximum(abs, u), maximum(abs, v), maximum(abs, w))

# === Timed measurement: hand off to the shared SYPD harness ===
measure_sypd!(simulation; name="eady_les",
              extra_config=(; Nx, Ny, Nz, Lx, Lz, weno_order,
                            initial_dt, max_dt, cfl, spinup_days,
                            spinup_wall_seconds=spinup_wall,
                            spinup_iters=spinup_iters,
                            spinup_final_dt=spinup_final_dt))
