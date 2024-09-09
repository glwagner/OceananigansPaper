#=
using Oceananigans

const δ = 10

# We start by introducing a Gaussian wave envelope with width `δ`,
A(ξ) = exp(- ξ^2 / 2δ^2)

# The envelope has derivatives
A′(ξ) = - ξ / δ^2 * A(ξ)
A′′(ξ) = ((ξ / δ)^2 - 1) * A(ξ) / δ^2

# Next, defining the Stokes drift as
uˢ(x, z, t) = A(x - t) * exp(z)

# we obtain the derivatives
∂z_uˢ(x, z, t) = + A(x - t) * exp(z)
∂t_uˢ(x, z, t) = - A′(x - t) * exp(z)

# Further using Stokes drift non-divergence (Vanneste and Young 2022), which
# requires ∂z_wˢ = - ∂x_uˢ, we find that wˢ = - A′(x - t) exp(z), and thus
∂x_wˢ(x, z, t) = - A′′(x - t) * exp(z)
∂t_wˢ(x, z, t) = + A′′(x - t) * exp(z)

stokes_drift = StokesDrift(; ∂z_uˢ, ∂t_uˢ, ∂t_wˢ, ∂x_wˢ)

grid = RectilinearGrid(size = (256, 128), x = (-100, 200), z = (-100, 0),
                       topology = (Periodic, Flat, Bounded))

model = NonhydrostaticModel(; grid, stokes_drift, timestepper = :RungeKutta3)

# Set Lagrangian-mean flow equal to uˢ,
uᵢ(x, z) = uˢ(x, z, 0)
set!(model, u=uᵢ)

Δx = minimum_xspacing(grid)
Δt = 0.2 * Δx
simulation = Simulation(model; Δt, stop_time = 100)

progress(sim) = @info string("Iter: ", iteration(sim), ", time: ", time(sim))
simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

filename = "surface_wave_induced_flow.jld2"
outputs = model.velocities
simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs; filename,
                                                    schedule = IterationInterval(10),
                                                    overwrite_existing = true)

run!(simulation)
=#

u, v, w = model.velocities
t = time(simulation)
x = xnodes(u)

using GLMakie

uˢ(x, z, t) = A(x - t) * exp(z)
uˢn = XFaceField(grid)
set!(uˢn, (x, z) -> uˢ(x, z, t))
u′ = Field(u - uˢn)
compute!(u′)

fig = Figure(size=(1200, 1200))
axA = Axis(fig[1, 1])
axw = Axis(fig[2, 1])
axu = Axis(fig[3, 1])
axu′ = Axis(fig[4, 1])

lines!(axA, x, A.(x .- t))
lines!(axA, x, A′.(x .- t))

heatmap!(axw, w)
#heatmap!(axu, u)
heatmap!(axu′, u′)

x, y, z = nodes(u)
Nx, Ny, Nz = size(u)

s = 1
arrows!(axu,
        x[1:s:Nx],
        z[1:s:Nz],
        interior(u, 1:s:Nx, 1, 1:s:Nz),
        interior(w, 1:s:Nx, 1, 1:s:Nz),
        lengthscale = 20)

for ax in (axw, axu, axu′)
    xlims!(axA, 50, 150)
    xlims!(ax, 50, 150)
    ylims!(ax, -20, 0)
end

#=
include("streamlines.jl")

x₀ = collect(70:10.0:130)
Np = length(x₀)

y₀ = zeros(Np)
z₀ = -0.1 * ones(Np)

initial_positions = (x₀, y₀, z₀)
samples = streamlines(model.velocities, initial_positions, 10000, 100)


α = 0.5
for p = 1:Np
    lines!(axu, samples[p].x, samples[p].z, color=(:black, α))
end

xlims!(axu, -100, 200)
ylims!(axu, -100, 0)
=#

display(fig)
