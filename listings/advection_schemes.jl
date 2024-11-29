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
    advect_tracer(tracer_advection, closure=nothing, Nx=128)

Advect a tracer with a top-hat profile until `stop_time` using `Nx` grid points.
"""
function advect_tracer(tracer_advection; closure=nothing, Nx=128, stop_time=1)
    grid = RectilinearGrid(size=Nx, x=(-1, 1), halo=6, topology=(Periodic, Flat, Flat))

    model = NonhydrostaticModel(; grid, closure, tracer_advection, tracers = :c)

    set!(model, u=1, c=cᵢ)
    simulation = Simulation(model; Δt=0.5/Nx, stop_time)
    run!(simulation)
    return model.tracers.c
end

using GLMakie

set_theme!(Theme(fontsize=18, linewidth=3, linealpha=0.6))

fig = Figure(size=(800, 240))
ax = Axis(fig[1, 1], xlabel="x") #, ylabel="c(t=4)")

stop_time = 0.1
c = advect_tracer(Centered(order=2); stop_time)
lines!(ax, c, label="t=0.1, Centered(order=2)")

c = advect_tracer(UpwindBiased(order=3); stop_time)
lines!(ax, c, label="t=0.1, UpwindBiased(order=3)")

#c = advect_tracer(UpwindBiased(order=5); stop_time)
#lines!(ax, c, label="UpwindBiased(order=5)")

stop_time = 4

for N = (5, 11)
    scheme = WENO(order=N)
    c = advect_tracer(scheme; stop_time)
    lines!(ax, c, label="t=4, WENO(order=$N)")
end

solution(x) = cᵢ(x)
set!(c, solution)
lines!(ax, c, linestyle=:dash, color=(:black, 0.2), label="Initial condition")

solution(x) = cᵢ(x-stop_time) > 1 ? 0 : 1 
set!(c, solution)
lines!(ax, c, linestyle=:dash, color=(:black, 0.6), label="Exact solution at t=4")

# c = advect_tracer(WENO(order=11), closure=ScalarDiffusivity(κ=0.04))
# lines!(ax, c, label="WENO(order=11) with κ=0.04")

xlims!(ax, -2, 6)
hidespines!(ax, :t, :r, :l)
hideydecorations!(ax, grid=false)

Legend(fig[1, 0], ax, framevisible=false)

display(fig)

save("advection_schemes.png", fig)

