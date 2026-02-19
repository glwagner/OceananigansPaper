using Oceananigans
using Oceananigans.Units
import NumericalEarth
using Printf

# Grid setup: 1/6° global resolution
Nz = 40
exponential_z_faces(; Nz, depth) = -depth * (1 .- (0:Nz) / Nz).^2
z = exponential_z_faces(; Nz, depth=5000)
arch = GPU()
grid = LatitudeLongitudeGrid(arch; z,
                             size = (2160, 1020, Nz),
                             halo = (7, 7, 7),
                             longitude = (0, 360),
                             latitude = (-70, 70))

bathymetry = NumericalEarth.regrid_bathymetry(grid)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry))

# Build coupled Earth System Model
ocean = NumericalEarth.ocean_simulation(grid)
set!(ocean.model, T=20, S=35)

atmosphere = NumericalEarth.JRA55PrescribedAtmosphere(arch)
coupled_model = NumericalEarth.OceanOnlyModel(ocean; atmosphere)
simulation = Simulation(coupled_model, Δt=10minutes, stop_iteration=100)

# Disable the NaN checker so benchmark runs to completion
pop!(simulation.callbacks, :nan_checker, nothing)

# Warmup: compile and warm up (no output writers)
@info "Warming up for 100 iterations..."
run!(simulation)
@info "Warmup complete."

# Helper to reset model state between trials
function reset_model!(simulation, ocean)
    simulation.model.clock.iteration = 0
    simulation.model.clock.time = 0
    set!(ocean.model, T=20, S=35)
    u, v, w = ocean.model.velocities
    set!(u, 0)
    set!(v, 0)
    set!(w, 0)
    return nothing
end

# Benchmark parameters
Nsteps = 1440
output_intervals = [
    ("Control (no I/O)", nothing),
    ("k=144 (daily)",    144),
    ("k=6 (hourly)",     6),
    ("k=2 (20 min)",     2),
    ("k=1 (every step)", 1),
]

# 3D ocean fields to write
benchmark_outputs = merge(ocean.model.velocities, ocean.model.tracers)

results = []

for (name, k) in output_intervals
    # Reset clock and fields
    reset_model!(simulation, ocean)
    simulation.stop_iteration = Nsteps

    # Attach output writer if this is not the control run
    benchmark_file = "io_benchmark_output.jld2"
    if !isnothing(k)
        writer = JLD2Writer(ocean.model, benchmark_outputs;
                            filename = benchmark_file,
                            array_type = Array{Float32},
                            schedule = IterationInterval(k),
                            overwrite_existing = true)
        simulation.output_writers[:benchmark] = writer
    end

    @info "Running trial: $name ($Nsteps iterations)..."
    elapsed = @elapsed run!(simulation)

    # Clean up output writer and file
    if !isnothing(k)
        pop!(simulation.output_writers, :benchmark)
        rm(benchmark_file, force=true)
    end

    push!(results, (; name, k, elapsed))
    @info @sprintf("  %s: %.2f s (%.4f s/iter)", name, elapsed, elapsed / Nsteps)
end

# Print results table
control_time = results[1].elapsed
simulated_days = Nsteps * 10 / (60 * 24)  # 10 min per step → days

println()
println("=" ^ 80)
println("I/O Benchmark Results — 1/6° global ESM, $Nsteps iterations ($(@sprintf("%.1f", simulated_days)) simulated days)")
println("=" ^ 80)
@printf("%-22s %10s %12s %8s %10s\n", "Trial", "Time (s)", "s/iteration", "SYPD", "Overhead")
println("-" ^ 80)

for r in results
    sypd = simulated_days / r.elapsed * 365
    overhead = (r.elapsed - control_time) / control_time * 100
    overhead_str = r.name == "Control (no I/O)" ? "—" : @sprintf("%.1f%%", overhead)
    @printf("%-22s %10.2f %12.4f %8.2f %10s\n",
            r.name, r.elapsed, r.elapsed / Nsteps, sypd, overhead_str)
end

println("=" ^ 80)
