using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Printf
using CUDA
using CairoMakie

function gaussian_islands_tripolar_grid(arch, Nx, Ny, Nz; halo=(7, 7, 7))
    underlying_grid = TripolarGrid(arch; size=(Nx, Ny, Nz), halo, z=(-3000, 0))

    zb = z[1]
    h = -zb + 100
    gaussian_islands(λ, φ) = zb + h * (mtn₁(λ, φ) + mtn₂(λ, φ))

    return ImmersedBoundaryGrid(underlying_grid,
                                GridFittedBottom(gaussian_islands);
                                active_cells_map = false)
end

grid_type = "lat_lon"
grid_type = "tripolar"
filename = "baroclinic_wave_" * grid_type

arch = GPU()

resolution = 1//2
Nx = 360 ÷ resolution
Ny = 160 ÷ resolution
Nz = 10

grid = if grid_type == "lat_lon"
    LatitudeLongitudeGrid(arch;
                          size = (Nx, Ny, Nz),
                          halo = (7, 7, 7),
                          latitude = (-80, 80),
                          longitude = (0, 360),
                          z = (-3000, 0))
elseif grid_type == "tripolar"
    gaussian_islands_tripolar_grid(arch, Nx, Ny, Nz)
end

momentum_advection = WENOVectorInvariant(order=9)
tracer_advection   = WENO(order=7)
coriolis = HydrostaticSphericalCoriolis()
buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())

model = HydrostaticFreeSurfaceModel(; grid, coriolis,
                                      buoyancy, tracers=(:T, :S),
                                      momentum_advection, tracer_advection)

Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 45) / 8)) / 2 + rand()
Sᵢ(λ, φ, z) = 28 - 5e-3 * z + rand()
set!(model, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=5minutes, stop_time=200days)

function progress(sim)
    T = sim.model.tracers.T
    S = sim.model.tracers.S
    u, v, w = sim.model.velocities

    msg = @sprintf("%d: %s, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim),
                   maximum(abs, u),
                   maximum(abs, v),
                   maximum(abs, w))

    msg *= @sprintf(", extrema(T): (%.2f, %.2f), extrema(S): (%.2f, %.2f)",
                    minimum(T), maximum(T), minimum(S), maximum(S))

    @info msg

    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
ζ = ∂x(v) - ∂y(u)
fields = (; u, v, w, T, S, ζ)

ow = JLD2OutputWriter(model, fields,
                      filename = filename * ".jld2",
                      indices = (:, :, grid.Nz),
                      schedule = TimeInterval(0.5days),
                      overwrite_existing = true)

simulation.output_writers[:surface] = ow

run!(simulation)

#=
using Oceananigans
using Oceananigans.Units
using GLMakie


u = FieldTimeSeries(filename * ".jld2", "u"; backend = OnDisk())
v = FieldTimeSeries(filename * ".jld2", "v"; backend = OnDisk())
T = FieldTimeSeries(filename * ".jld2", "T"; backend = OnDisk())
ζ = FieldTimeSeries(filename * ".jld2", "ζ"; backend = OnDisk())

times = u.times
Nt = length(times)

n = Observable(Nt)

ζn = @lift ζ[$n]

Tn = @lift T[$n]
Sn = @lift S[$n]

un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)
s = Field(sqrt(un^2 + vn^2))

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
end

fig = Figure(size=(800, 1200))

kwargs_axis = (xticks = 0:60:360, yticks = -80:40:80)
axs = Axis(fig[2, 1]; kwargs_axis...)
axζ = Axis(fig[3, 1]; kwargs_axis...)
axT = Axis(fig[4, 1]; kwargs_axis...)

hm = heatmap!(axs, sn, colorrange=(0, 6))
Colorbar(fig[2, 2], hm, label = "Surface speed (m s⁻¹)", labelsize=20)

hm = heatmap!(axζ, ζn, colormap=:balance, colorrange=(-3e-5, 3e-5))
ticks = (-2e-5:1e-5:2e-5, ["-2", "-1", "0", "1", "2"])
Colorbar(fig[3, 2], hm; ticks, label = "relative vorticity (10⁻⁵ s⁻¹)", labelsize=20)

hm = heatmap!(axT, Tn, colormap=:thermal, colorrange=(0, 31))
Colorbar(fig[4, 2], hm, label = "surface temperature (ᵒC)", labelsize=20)

title = @lift @sprintf("%s days", round(times[$n] / day, digits=2))
fig[1, :] = Label(fig, title, tellwidth=false)

fig


n[] = 1
save(filename * "_initial_snapshot.png", fig)

n[] = Nt
save(filename * "_final_snapshot.png", fig)

record(fig, filename * ".mp4", 1:Nt, framerate = 16) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end


function geographic2cartesian(λ, φ, r=1)
    Nλ = length(λ)
    Nφ = length(φ)
    λ = repeat(reshape(λ, Nλ, 1), 1, Nφ)
    φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end

u = FieldTimeSeries(filename * ".jld2", "u"; backend = OnDisk())
v = FieldTimeSeries(filename * ".jld2", "v"; backend = OnDisk())
T = FieldTimeSeries(filename * ".jld2", "T"; backend = OnDisk())
ζ = FieldTimeSeries(filename * ".jld2", "ζ"; backend = OnDisk())

times = u.times
Nt = length(times)

grid = u.grid
λ = λnodes(grid, Center(), Center(), Center())
φ = φnodes(grid, CeP0+r\nter(), Center(), Center())
x, y, z = geographic2cartesian(λ, φ)

Nx, Ny, Nz = size(grid)

n = Observable(1)

ζn = @lift 1e5 * interior(ζ[$n], :, :, 1)
Tn = @lift interior(T[$n], :, :, 1)

un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)
s = Field(sqrt(un^2 + vn^2))

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
    interior(s, :, :, 1)
end

fig = Figure(size=(1400, 700))

kw = (elevation= deg2rad(-10), azimuth=deg2rad(90), aspect=:equal)

axs = Axis3(fig[1, 1]; kw...)
axζ = Axis3(fig[1, 2]; kw...)
axT = Axis3(fig[1, 3]; kw...)

sf = surface!(axs, x, y, z, color=sn, colorrange=(0, 6), colormap=:viridis, nan_color=:lightgray)
Colorbar(fig[2, 1], sf, width=Relative(0.6), vertical=false, label = "Surface speed (m s⁻¹)", labelsize=20)

sf = surface!(axζ, x, y, z, color=ζn, colorrange=(-3, 3), colormap=:balance, nan_color=:lightgray)
Colorbar(fig[2, 2], sf, width=Relative(0.6), vertical=false, label = "relative vorticity (10⁻⁵ s⁻¹)", labelsize=20)

sf = surface!(axT, x, y, z, color=Tn, colorrange=(0, 31), colormap=:thermal, nan_color=:lightgray)
Colorbar(fig[2, 3], sf, width=Relative(0.6), vertical=false, label = "surface temperature (ᵒC)", labelsize=20)


t = u.times
title = @lift @sprintf("%s days", round((times[$n]+ 1e-6) / day, digits=2))
Label(fig[0, :], title, fontsize=24)

hidedecorations!(axζ)
hidedecorations!(axs)
hidedecorations!(axT)
hidespines!(axζ)
hidespines!(axs)
hidespines!(axT)

rowgap!(fig.layout, 1, Relative(-0.1))
colgap!(fig.layout, 1, Relative(-0.07))
colgap!(fig.layout, 2, Relative(-0.07))

fig

n[] = 1
save(filename * "_initial_snapshot_sphere.png", fig)

n[] = Nt
save(filename * "_final_snapshot_sphere.png", fig)

n[] = 1
save(filename * "_final_snapshot_sphere.png", fig)

record(fig, filename * "_sphere.mp4", 1:Nt, framerate = 16) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
=#
