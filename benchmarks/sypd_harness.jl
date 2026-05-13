# Shared SYPD measurement harness.
# A listing or benchmark script builds a `simulation` object, then calls
# `measure_sypd!(simulation; warmup, timed, name)` to get a stable SYPD reading.
# The harness:
#   1. Clears output writers + diagnostics so I/O isn't on the timing path.
#   2. Runs `warmup` iterations to flush JIT / let CFL settle.
#   3. Times the next `timed` iterations.
#   4. Computes SYPD = (timed × mean_Δt) / (365.25 × wall_seconds).
#   5. Writes a JSON summary to benchmarks/runs/<name>_<jobid>.json.
#
# Override knobs (env vars):
#   BENCHMARK_WARMUP_ITERS   — default 50
#   BENCHMARK_TIMED_ITERS    — default 200

using Printf
using Dates

function _write_json(io, x; indent=0)
    pad = "  " ^ indent
    if x isa NamedTuple
        println(io, "{")
        ks = collect(keys(x))
        for (i, k) in enumerate(ks)
            print(io, "  " ^ (indent+1), "\"", k, "\": ")
            _write_json(io, x[k]; indent=indent+1)
            i < length(ks) ? println(io, ",") : println(io)
        end
        print(io, pad, "}")
    elseif x isa AbstractVector
        if isempty(x); print(io, "[]")
        else
            println(io, "[")
            for (i, v) in enumerate(x)
                print(io, "  " ^ (indent+1)); _write_json(io, v; indent=indent+1)
                i < length(x) ? println(io, ",") : println(io)
            end
            print(io, pad, "]")
        end
    elseif x isa AbstractString
        print(io, "\"", replace(x, "\\" => "\\\\", "\"" => "\\\""), "\"")
    elseif x isa Bool; print(io, x ? "true" : "false")
    elseif x isa Number; print(io, x)
    else; print(io, "\"", string(x), "\"")
    end
end

function measure_sypd!(simulation; warmup::Int=parse(Int, get(ENV, "BENCHMARK_WARMUP_ITERS", "50")),
                                    timed::Int=parse(Int, get(ENV, "BENCHMARK_TIMED_ITERS", "200")),
                                    name::String="benchmark",
                                    extra_config::NamedTuple=NamedTuple())

    empty!(simulation.output_writers)
    empty!(simulation.diagnostics)

    sim_t0_warmup = simulation.model.clock.time
    iter_t0_warmup = simulation.model.clock.iteration
    @info @sprintf("[%s] Warmup: %d iterations starting at iter=%d, sim_t=%.3g s",
                   name, warmup, iter_t0_warmup, sim_t0_warmup)
    simulation.stop_iteration = iter_t0_warmup + warmup
    simulation.stop_time = Inf
    wall_start_warmup = time()
    run!(simulation)
    wall_warmup = time() - wall_start_warmup
    sim_after_warmup = simulation.model.clock.time
    iter_after_warmup = simulation.model.clock.iteration
    @info @sprintf("[%s] Warmup done: %d iters in %.2f s, Δt~%.4g s",
                   name, iter_after_warmup - iter_t0_warmup, wall_warmup,
                   (sim_after_warmup - sim_t0_warmup) / max(1, iter_after_warmup - iter_t0_warmup))

    @info @sprintf("[%s] Timed run: %d iterations starting at iter=%d, sim_t=%.3g s",
                   name, timed, iter_after_warmup, sim_after_warmup)
    simulation.stop_iteration = iter_after_warmup + timed
    wall_start_timed = time()
    run!(simulation)
    wall_timed = time() - wall_start_timed
    sim_after_timed = simulation.model.clock.time
    iter_after_timed = simulation.model.clock.iteration

    iters_done   = iter_after_timed - iter_after_warmup
    sim_seconds  = sim_after_timed  - sim_after_warmup
    mean_dt      = sim_seconds / max(1, iters_done)
    sec_per_iter = wall_timed / max(1, iters_done)
    sypd         = sim_seconds / (365.25 * 86400) / (wall_timed / 86400)
    sdph         = (sim_seconds / 86400) / (wall_timed / 3600)

    @info @sprintf("[%s] Timed done: %d iters in %.2f s wall, mean_Δt=%.4g s",
                   name, iters_done, wall_timed, mean_dt)
    @info @sprintf("[%s] sec/iter=%.4f  SDPH=%.3f sim_days/hr  SYPD=%.4f",
                   name, sec_per_iter, sdph, sypd)

    gpu_name = try
        # Caller is responsible for `using CUDA` before include()ing the harness.
        String(CUDA.name(CUDA.device()))
    catch; "unknown"; end

    result = (
        timestamp = string(now()),
        name = name,
        config = extra_config,
        hardware = (
            gpu = gpu_name,
            n_gpus = 1,
            n_nodes = 1,
            slurm_jobid = get(ENV, "SLURM_JOB_ID", "n/a"),
            slurm_nodelist = get(ENV, "SLURM_JOB_NODELIST", "n/a"),
        ),
        warmup = (
            iterations = iter_after_warmup - iter_t0_warmup,
            wall_seconds = wall_warmup,
        ),
        timed = (
            iterations = iters_done,
            wall_seconds = wall_timed,
            mean_dt_seconds = mean_dt,
            sec_per_iter = sec_per_iter,
            sim_seconds = sim_seconds,
            sim_days = sim_seconds / 86400,
            sdph = sdph,
            sypd = sypd,
        ),
        versions = (
            julia = string(VERSION),
            oceananigans = string(pkgversion(parentmodule(typeof(simulation.model)))),
        ),
    )

    runs_dir = joinpath(@__DIR__, "runs")
    mkpath(runs_dir)
    outpath = joinpath(runs_dir, "$(name)_$(get(ENV, "SLURM_JOB_ID", "local")).json")
    open(outpath, "w") do io
        _write_json(io, result); println(io)
    end
    @info "[$name] wrote $outpath"
    return result
end
