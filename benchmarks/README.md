# Benchmarks

Short SYPD (simulated years per day) measurements for the listings the JAMES
revision needs timing for. Each benchmark builds the model at the paper-stated
configuration, warms the GPU for `BENCHMARK_WARMUP_ITERS` iterations (default
50) to flush JIT and let CFL settle, then times the next `BENCHMARK_TIMED_ITERS`
iterations (default 200) and computes:

```
SYPD = (timed × mean_Δt) / (365.25 × wall_seconds)
```

Output writers and diagnostics are cleared before timing so I/O isn't on the
hot path.

Each run produces a JSON summary at `benchmarks/runs/<name>_<jobid>.json`.

## Files

- `sypd_harness.jl` — shared helper. Defines `measure_sypd!(simulation; ...)`.
- `eady_les_benchmark.jl` — Sim 1 (paper config: 1000×1000×64, WENO(9), single H100).
- `flow_past_headland_benchmark.jl` — Sim 3 (paper config: 384×128×64, WENO(9), single H100).
- `global_simulation_benchmark.jl` — global ocean–sea-ice at the listing's 1° config (1440×560×40, single H100). Needs `ECCO_USERNAME` + `ECCO_WEBDAV_PASSWORD` env vars.
- `runs/` — per-job JSON summaries (gitignored).

## How to run

From the repo root, on the AWS pcluster head node:

```bash
# Eady
BENCH_FILE=benchmarks/eady_les_benchmark.jl \
  sbatch --job-name=bench-eady \
         --export=ALL,BENCH_FILE=benchmarks/eady_les_benchmark.jl \
         scripts/run_benchmark.batch

# Headland
BENCH_FILE=benchmarks/flow_past_headland_benchmark.jl \
  sbatch --job-name=bench-headland \
         --export=ALL,BENCH_FILE=benchmarks/flow_past_headland_benchmark.jl \
         scripts/run_benchmark.batch

# Global (needs ECCO creds in ~/.secrets/ecco.env; the batch script sources them)
BENCH_FILE=benchmarks/global_simulation_benchmark.jl \
  sbatch --job-name=bench-global \
         --export=ALL,BENCH_FILE=benchmarks/global_simulation_benchmark.jl \
         scripts/run_benchmark.batch
```

`scripts/run_benchmark.batch` pins `--nodelist=gpu-prod-st-gpu-prod-2` and uses
a single H100 (`--gres=gpu:1`).

## Tunables (env vars)

| Variable | Default | Effect |
|---|---|---|
| `BENCHMARK_WARMUP_ITERS` | 50 | iterations run untimed before measurement |
| `BENCHMARK_TIMED_ITERS` | 200 | iterations measured |
| `EADY_NX`, `EADY_NY`, `EADY_NZ` | 1000, 1000, 64 | Eady grid |
| `EADY_LX`, `EADY_LZ` | 4000, 128 (m) | Eady domain |
| `EADY_WENO` | 9 | Eady advection order |
| `HEADLAND_NZ` | 64 | Headland vertical resolution (grid scales 6Nz × 2Nz × Nz) |
| `GLOBAL_NX`, `GLOBAL_NY`, `GLOBAL_NZ` | 1440, 560, 40 | global grid |
| `GLOBAL_DT_MIN` | 5 | global Δt in minutes |

## Output JSON shape

```json
{
  "name": "eady_les",
  "config": { "Nx": 1000, "Ny": 1000, "Nz": 64, ... },
  "hardware": { "gpu": "NVIDIA H100 80GB HBM3", "n_gpus": 1, "n_nodes": 1, ... },
  "warmup": { "iterations": 50, "wall_seconds": 12.3 },
  "timed":  { "iterations": 200, "wall_seconds": 41.2,
              "mean_dt_seconds": 28.7, "sec_per_iter": 0.206,
              "sim_seconds": 5740, "sim_days": 0.066,
              "sdph": 5.78, "sypd": 0.01585 },
  "versions": { "julia": "1.12.6", "oceananigans": "0.107.4" }
}
```

## Reproducibility

The benchmark scripts pin `JULIA_CPU_TARGET=generic` so the precompile cache is
portable across nodes. The Manifest.toml on the branch this work landed on
(`glw/update-packages-2026-05`) records the exact package versions used.
