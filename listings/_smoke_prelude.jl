"""
Opt-in smoke-test helper for listings.

When `OCEANANIGANS_PAPER_SMOKE=1`, `smoke_test_simulation!(simulation)` clamps
the simulation to a small number of iterations and clears output writers so
the listing runs end-to-end in seconds. When the env var is unset, this is a
no-op and the listing behaves exactly as in the paper.

Optional: `OCEANANIGANS_PAPER_SMOKE_ITERS` (default 10) sets the iteration cap.
"""
function smoke_test_simulation!(simulation; iters::Int=10)
    if get(ENV, "OCEANANIGANS_PAPER_SMOKE", "0") == "1"
        n = parse(Int, get(ENV, "OCEANANIGANS_PAPER_SMOKE_ITERS", string(iters)))
        target = simulation.model.clock.iteration + n
        simulation.stop_iteration = min(simulation.stop_iteration, target)
        empty!(simulation.output_writers)
        empty!(simulation.diagnostics)
        @info "[SMOKE] stop_iteration clamped to $(simulation.stop_iteration); output writers cleared"
    end
    return simulation
end
