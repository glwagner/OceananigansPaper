#=
using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Printf
using CUDA

# grid_type = "lat_lon"
# grid_type = "tripolar"
grid_type = "cubed_sphere"
filename = "baroclinic_wave_" * grid_type

arch = CPU()
res = 1//2 # degrees
Nx, Ny, Nz = 360 ÷ res, 170 ÷ res, 10
H = 3000 # domain depth

# Simple lat-lon grid
grid = if grid_type == "lat_lon"
    LatitudeLongitudeGrid(arch; size=(Nx, Ny, Nz), halo=(7, 7, 7),
                          latitude=(-85, 85), longitude=(0, 360), z=(-H, 0))
elseif grid_type == "tripolar"
    underlying_grid = TripolarGrid(arch; size=(Nx, Ny, Nz), halo=(7, 7, 7), z=(-H, 0))

    # Add bathymetry to cover north pole singularities:
    dφ, dλ, λ₀, φ₀, h = 10, 20, 70, 55, 100
    isle(λ, φ) = exp(-(λ - λ₀)^2 / 2dλ^2 - (φ - φ₀)^2 / 2dφ^2)
    gaussian_isles(λ, φ) = - H + (H + h) * (isle(λ, φ) + isle(λ - 180, φ))

    ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(gaussian_isles))
elseif grid_type == "cubed_sphere"
    Nh = 32 #Int(Nx / 6)
    ConformalCubedSphereGrid(arch, panel_size=(Nh, Nh, Nz), horizontal_direction_halo=7, z=(-H, 0))
else
    error("Not sure what grid we should make for grid_type=$grid_type !")
end
 
momentum_advection = WENOVectorInvariant(order=9)
tracer_advection = WENO(order=7)
coriolis = HydrostaticSphericalCoriolis()
buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())
free_surface = SplitExplicitFreeSurface(grid, substeps=30)

model = HydrostaticFreeSurfaceModel(; grid, coriolis, free_surface,
                                      buoyancy, tracers=(:T, :S),
                                      momentum_advection, tracer_advection)

Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 45) / 8)) / 2 + rand()
Sᵢ(λ, φ, z) = 28 - 5e-3 * z + rand()
set!(model, T=Tᵢ, S=Sᵢ)

#simulation = Simulation(model, Δt=5minutes, stop_time=200days)
simulation = Simulation(model, Δt=5minutes, stop_iteration=10)

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

add_callback!(simulation, progress, IterationInterval(1))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
ζ = ∂x(v) - ∂y(u)
fields = (; u, v, w, T, S, ζ)

ow = JLD2Writer(model, fields,
                filename = filename * ".jld2",
                indices = (:, :, grid.Nz),
                schedule = TimeInterval(0.5days),
                overwrite_existing = true)

simulation.output_writers[:surface] = ow

run!(simulation)
=#

using Oceananigans
using Oceananigans.Units
using CUDA
using GLMakie
using Printf

function geographic2cartesian(λ, φ, r=1)
    if ndims(λ) == 1
        Nλ = length(λ)
        Nφ = length(φ)
        λ = repeat(reshape(λ, Nλ, 1), 1, Nφ)
        φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)
    end

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end

fig = Figure(size=(1200, 900))

kw = (elevation= deg2rad(50), azimuth=deg2rad(0), aspect=:equal)

axζ1 = Axis3(fig[1, 1]; kw...)
axT1 = Axis3(fig[2, 1]; kw...)
axζ2 = Axis3(fig[1, 2]; kw...)
axT2 = Axis3(fig[2, 2]; kw...)
axζ3 = Axis3(fig[1, 3]; kw...)
axT3 = Axis3(fig[2, 3]; kw...)

#slider = Slider(fig[3, 1:3], range=1:Nt, startvalue=1)
#n = slider.value
n = Observable(361)

filename1 = "baroclinic_wave_lat_lon"
filename2 = "baroclinic_wave_tripolar"
#filename3 = "baroclinic_wave_tripolar"
filename3 = "baroclinic_wave_rotated_lat_lon"

T1 = FieldTimeSeries(filename1 * ".jld2", "T"; backend = OnDisk())
ζ1 = FieldTimeSeries(filename1 * ".jld2", "ζ"; backend = OnDisk())

T2 = FieldTimeSeries(filename2 * ".jld2", "T"; backend = OnDisk())
ζ2 = FieldTimeSeries(filename2 * ".jld2", "ζ"; backend = OnDisk())

T3 = FieldTimeSeries(filename3 * ".jld2", "T") #; backend = OnDisk())
ζ3 = FieldTimeSeries(filename3 * ".jld2", "ζ") #; backend = OnDisk())

times = T1.times
Nt = length(times)

grid1 = T1.grid
λ1 = λnodes(grid1, Center(), Center(), Center())
φ1 = φnodes(grid1, Center(), Center(), Center())
x1, y1, z1 = geographic2cartesian(λ1, φ1)

grid2 = T2.grid
λ2 = λnodes(grid2, Center(), Center(), Center())
φ2 = φnodes(grid2, Center(), Center(), Center())
x2, y2, z2 = geographic2cartesian(λ2, φ2)

# Special treatment for cubed sphere stuff
grid3 = LatitudeLongitudeGrid(size=(4 * 360, 4 * 180, 1), longitude=(0, 360), latitude=(-90, 90), z=(-1, 0))
# grid3 = T3.grid
λ3 = λnodes(grid3, Center(), Center(), Center())
φ3 = φnodes(grid3, Center(), Center(), Center())
x3, y3, z3 = geographic2cartesian(λ3, φ3)

ζ1n = @lift 1e5 * interior(ζ1[$n], :, :, 1)
T1n = @lift interior(T1[$n], :, :, 1)

ζ2n = @lift 1e5 * interior(ζ2[$n], :, :, 1)
T2n = @lift interior(T2[$n], :, :, 1)

ζ3llg = Field{Face, Face, Center}(grid3)
T3llg = CenterField(grid3)

using Oceananigans.Fields: interpolate!

ζ3n = @lift begin
    interpolate!(ζ3llg, ζ3[$n])
    1e5 * interior(ζ3llg, :, :, 1)
end

T3n = @lift begin
    interpolate!(T3llg, T3[$n])
    interior(T3llg, :, :, 1)
end

sf = surface!(axζ1, x1, y1, z1, color=ζ1n, colorrange=(-3, 3), colormap=:balance, nan_color=:lightgray)
sf = surface!(axζ2, x2, y2, z2, color=ζ2n, colorrange=(-3, 3), colormap=:balance, nan_color=:lightgray)
sf = surface!(axζ3, x3, y3, z3, color=ζ3n, colorrange=(-3, 3), colormap=:balance, nan_color=:lightgray)

Colorbar(fig[1, 4], sf, height=Relative(0.6), flipaxis=true, vertical=true,
         label="Relative vorticity (10⁻⁵ s⁻¹)", labelsize=20)

sf = surface!(axT1, x1, y1, z1, color=T1n, colorrange=(0, 31), colormap=:thermal, nan_color=:lightgray)
sf = surface!(axT2, x2, y2, z2, color=T2n, colorrange=(0, 31), colormap=:thermal, nan_color=:lightgray)
sf = surface!(axT3, x3, y3, z3, color=T3n, colorrange=(0, 31), colormap=:thermal, nan_color=:lightgray)

Colorbar(fig[2, 4], sf, height=Relative(0.6), flipaxis=true, vertical=true,
         label="Surface temperature (ᵒC)", labelsize=20)

title = @lift @sprintf("%s days", round((times[$n]+ 1e-6) / day, digits=2))
#Label(fig[0, :], title, fontsize=24)

Label(fig[0, 1], "Latitude-Longitude Grid", tellwidth=false, fontsize=24)
Label(fig[0, 2], "Tripolar Grid", tellwidth=false, fontsize=24)
Label(fig[0, 3], "Cubed Sphere Grid", tellwidth=false, fontsize=24)

hidedecorations!(axζ1)
hidedecorations!(axT1)
hidedecorations!(axζ2)
hidedecorations!(axT2)
hidedecorations!(axζ3)
hidedecorations!(axT3)
hidespines!(axζ1)
hidespines!(axT1)
hidespines!(axζ2)
hidespines!(axT2)
hidespines!(axζ3)
hidespines!(axT3)

rowgap!(fig.layout, 1, Relative(-0.1))
rowgap!(fig.layout, 2, Relative(-0.2))
colgap!(fig.layout, 1, Relative(-0.05))
colgap!(fig.layout, 2, Relative(-0.05))

fig

#=
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
