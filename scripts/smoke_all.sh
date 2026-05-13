#!/bin/bash
# Submit a smoke job for every runnable listing.
# Usage: scripts/smoke_all.sh [--gpu-only|--cpu-only|--include glob]
set -euo pipefail
cd "$(dirname "$0")/.."

# Listings that build a GPU model (need a GPU node).
GPU_LISTINGS=(
  baroclinic_wave
  cabbeling
  eady_les
  flow_around_cylinder
  flow_past_headland
  simple_global_simulation
)

# CPU-only listings without GLMakie that we can smoke now.
CPU_LISTINGS=(
  circling_tracer
  random_topography_internal_tide
  simple_advection
  two_dimensional_turbulence
)

# Illustrative-only listings (paper "fragments" with undefined refs to ax/fig/etc).
# Skipped from smoke; if you regenerate the paper figure, run them inside the
# corresponding figure script's notebook context.
NORUN_LISTINGS=()

# Listings deferred until we can run GLMakie on a headless node.
GLMAKIE_LISTINGS=(
  advection_schemes
  finite_volumes
  vertical_mixing_parameterizations
  vertically_implicit_diffusion
)

mode="${1:-all}"

submit_gpu() {
  local name="$1"
  echo "GPU  smoke-$name"
  LISTING="listings/${name}.jl" sbatch \
    --job-name="smoke-$name" \
    --export=ALL,LISTING="listings/${name}.jl" \
    scripts/smoke_one.batch
}

submit_cpu() {
  local name="$1"
  echo "CPU  smoke-$name"
  LISTING="listings/${name}.jl" sbatch \
    --job-name="smoke-$name" \
    --export=ALL,LISTING="listings/${name}.jl" \
    scripts/smoke_cpu.batch
}

case "$mode" in
  --gpu-only)
    for n in "${GPU_LISTINGS[@]}"; do submit_gpu "$n"; done
    ;;
  --cpu-only)
    for n in "${CPU_LISTINGS[@]}" "${NORUN_LISTINGS[@]}"; do submit_cpu "$n"; done
    ;;
  all|"")
    for n in "${GPU_LISTINGS[@]}"; do submit_gpu "$n"; done
    for n in "${CPU_LISTINGS[@]}" "${NORUN_LISTINGS[@]}"; do submit_cpu "$n"; done
    ;;
  *)
    echo "usage: $0 [--gpu-only|--cpu-only|all]" >&2
    exit 1
    ;;
esac

echo
echo "Deferred (need GLMakie on a display): ${GLMAKIE_LISTINGS[*]}"
echo
squeue -u "$USER" -o '%.10i %.20j %.10P %.10T %.10M %.6D %R'
