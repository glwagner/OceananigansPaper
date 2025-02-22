using Oceananigans

# the initial condition
@inline G(x, β, z) = exp(-β*(x - z)^2)
@inline F(x, α, a) = √(max(1 - α^2*(x-a)^2, 0.0))

const Z = -0.7
const δ = 0.005
const β = log(2)/(36*δ^2)
const a = 0.5
const α = 10

@inline function cᵢ(x)
    if x <= -0.6 && x >= -0.8
        return 1/6*(G(x, β, Z-δ) + 4*G(x, β, Z) + G(x, β, Z+δ))
    elseif x <= -0.2 && x >= -0.4
        return 1.0
    elseif x <= 0.2 && x >= 0.0
        return 1.0 - abs(10 * (x - 0.1))
    elseif x <= 0.6 && x >= 0.4
        return 1/6*(F(x, α, a-δ) + 4*F(x, α, a) + F(x, α, a+δ))
    else
        return 0.0
    end
end

"""
    advect_tracer(advection, closure=nothing, Nx=128)

Advect a tracer profile until `stop_time` using `Nx` grid points.
"""
function advect_tracer(advection; closure=nothing, Nx=128, stop_time=1)
    grid = RectilinearGrid(size=Nx, x=(-1, 1), halo=6, topology=(Periodic, Flat, Flat))

    model = NonhydrostaticModel(; grid, closure, advection, tracers = :c)

    set!(model, u=1, c=cᵢ)
    simulation = Simulation(model; Δt=0.7/Nx, stop_time)
    run!(simulation)
    return model.tracers.c
end

using GLMakie

set_theme!(Theme(fontsize=18, linewidth=3, linealpha=0.6))

fig = Figure(size=(1400, 500))
ax1 = Axis(fig[1, 1], xlabel="x", title="Linear reconstruction", yticks=[-1.5, -1, 0, 1, 1.5])  

stop_time = 2.0
c = advect_tracer(Centered(order=2); stop_time)
lines!(ax1, c, label="t=$(stop_time), Centered(order=2)", color = :palegreen3)

# If we want to show that this is equivalent to UpwindBiased(order=1)
# c = advect_tracer(Centered(order=2); closure = ScalarDiffusivity(κ=2/256), stop_time)
# scatter!(ax1, c, label="t=$(stop_time), Centered(order=2) + diffusion", color = :palegreen3)

c = advect_tracer(UpwindBiased(order=1); stop_time)
lines!(ax1, c, label="t=$(stop_time), UpwindBiased(order=1)", color = :salmon1)

c = advect_tracer(UpwindBiased(order=3); stop_time)
lines!(ax1, c, label="t=$(stop_time), UpwindBiased(order=3)", color = :firebrick1)

solution(x) = cᵢ(x)
set!(c, solution)
lines!(ax1, c, linestyle=:dash, color=(:black, 0.2))
xlims!(ax1, -1, 1)
ylims!(ax1, -0.5, 1.5)

ax2 = Axis(fig[1, 2], xlabel="x", title="WENO reconstruction", yticks=([-1.5, -1, 0, 1, 1.5], ["", "", "", "", ""]))
stop_time = 2.0

scheme = WENO(order=3)
c = advect_tracer(scheme; stop_time)
lines!(ax2, c, label="t=$(stop_time), WENO(order=3)", color=:dodgerblue3)

scheme = WENO(order=9)
c = advect_tracer(scheme; stop_time)
lines!(ax2, c, label="t=$(stop_time), WENO(order=9)", color=:lightskyblue)

solution(x) = cᵢ(x)
set!(c, solution)
lines!(ax2, c, linestyle=:dash, color=(:black, 0.2), label="Exact solution at t=$(stop_time)")

xlims!(ax2, -1,  1)
ylims!(ax2, -0.5, 1.5)

plots_in_fig  = AbstractPlot[]
labels_in_fig = AbstractString[]
for ax in (ax1, ax2)
    pl, lb = Makie.get_labeled_plots(ax, merge=false, unique=false)
    append!(plots_in_fig, pl)
    append!(labels_in_fig, lb)
end

ulabels = Base.unique(labels_in_fig)
mergedplots = [[lp for (i, lp) in enumerate(plots_in_fig) if labels_in_fig[i] == ul] for ul in ulabels]

Legend(fig[1, 0], mergedplots, ulabels)

display(fig)

save("advection_schemes.png", fig)

