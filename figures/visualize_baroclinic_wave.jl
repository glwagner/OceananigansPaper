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

kw = (elevation=deg2rad(70), azimuth=deg2rad(180), aspect=:equal)

axζ1 = Axis3(fig[1, 1]; kw...)
axT1 = Axis3(fig[2, 1]; kw...)
axζ2 = Axis3(fig[1, 2]; kw...)
axT2 = Axis3(fig[2, 2]; kw...)
axζ3 = Axis3(fig[1, 3]; kw...)
axT3 = Axis3(fig[2, 3]; kw...)

#slider = Slider(fig[3, 1:3], range=1:Nt, startvalue=1)
#n = slider.value
n = Observable(361)

filename1 = "baroclinic_wave_lat_lon_4"
filename2 = "baroclinic_wave_tripolar_4"
filename3 = "baroclinic_wave_rotated_lat_lon_4"

T1 = FieldTimeSeries(filename1 * ".jld2", "T"; backend = OnDisk())
ζ1 = FieldTimeSeries(filename1 * ".jld2", "ζ"; backend = OnDisk())
∇T1 = FieldTimeSeries(filename1 * ".jld2", "∇T"; backend = OnDisk())

T2 = FieldTimeSeries(filename2 * ".jld2", "T"; backend = OnDisk())
ζ2 = FieldTimeSeries(filename2 * ".jld2", "ζ"; backend = OnDisk())
∇T2 = FieldTimeSeries(filename2 * ".jld2", "∇T"; backend = OnDisk())

T3 = FieldTimeSeries(filename3 * ".jld2", "T") #; backend = OnDisk())
ζ3 = FieldTimeSeries(filename3 * ".jld2", "ζ") #; backend = OnDisk())
∇T3 = FieldTimeSeries(filename3 * ".jld2", "∇T"; backend = OnDisk())

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
# grid3 = LatitudeLongitudeGrid(size=(4 * 360, 4 * 180, 1), longitude=(0, 360), latitude=(-90, 90), z=(-1, 0))
grid3 = T3.grid
λ3 = λnodes(grid3, Center(), Center(), Center())
φ3 = φnodes(grid3, Center(), Center(), Center())
x3, y3, z3 = geographic2cartesian(λ3, φ3)

ζ1n = @lift 1e5 * interior(ζ1[$n], :, :, 1)
T1n = @lift interior(T1[$n], :, :, 1)
∇T1n = @lift interior(∇T1[$n], :, :, 1)

ζ2n = @lift begin
    ζ2i = 1e5 * interior(ζ2[$n], :, :, 1)
    ζ2i[ζ2i .== 0] .= NaN
    ζ2i
end

T2n = @lift begin
    T2i = interior(T2[$n], :, :, 1)
    T2i[T2i .== 0] .= NaN
    T2i
end

∇T2n = @lift interior(∇T2[$n], :, :, 1)

ζ3n = @lift 1e5 * interior(ζ3[$n], :, :, 1)
T3n = @lift interior(T3[$n], :, :, 1)
∇T3n = @lift interior(∇T3[$n], :, :, 1)

#=
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
=#

global_grid = LatitudeLongitudeGrid(size=(360, 180, 1), longitude=(0, 360), latitude=(-90, 90), z=(0, 1))
λg = λnodes(global_grid, Center(), Center(), Center())
φg = φnodes(global_grid, Center(), Center(), Center())
xg, yg, zg = geographic2cartesian(λg, φg, 0.95)
bg_field = CenterField(global_grid)
set!(bg_field, NaN)
bg = interior(bg_field, :, :, 1)

nan_color = :black
surface!(axζ1, xg, yg, zg, color=bg, colorrange=(-4, 4), colormap=:balance; nan_color)
surface!(axζ2, xg, yg, zg, color=bg, colorrange=(-4, 4), colormap=:balance; nan_color)
surface!(axζ3, xg, yg, zg, color=bg, colorrange=(-4, 4), colormap=:balance; nan_color)
sf = surface!(axζ1, x1, y1, z1, color=ζ1n, colorrange=(-4, 4), colormap=:balance; nan_color)
sf = surface!(axζ2, x2, y2, z2, color=ζ2n, colorrange=(-4, 4), colormap=:balance; nan_color)
sf = surface!(axζ3, x3, y3, z3, color=ζ3n, colorrange=(-4, 4), colormap=:balance; nan_color)

Colorbar(fig[1, 4], sf, height=Relative(0.6), flipaxis=true, vertical=true,
         label="Relative vorticity (10⁻⁵ s⁻¹)", labelsize=20)

surface!(axT1, xg, yg, zg, color=bg, colorrange=(2, 30), colormap=:thermal; nan_color)
surface!(axT2, xg, yg, zg, color=bg, colorrange=(2, 30), colormap=:thermal; nan_color)
surface!(axT3, xg, yg, zg, color=bg, colorrange=(2, 30), colormap=:thermal; nan_color)
sf = surface!(axT1, x1, y1, z1, color=T1n, colorrange=(2, 30), colormap=:thermal; nan_color)
sf = surface!(axT2, x2, y2, z2, color=T2n, colorrange=(2, 30), colormap=:thermal; nan_color)
sf = surface!(axT3, x3, y3, z3, color=T3n, colorrange=(2, 30), colormap=:thermal; nan_color)

Colorbar(fig[2, 4], sf, height=Relative(0.6), flipaxis=true, vertical=true,
         label="Surface temperature (ᵒC)", labelsize=20)

# sf = surface!(axT1, x1, y1, z1, color=∇T1n, colorrange=(0, 1e-4), colormap=:magma, nan_color=:lightgray)
# sf = surface!(axT2, x2, y2, z2, color=∇T2n, colorrange=(0, 1e-4), colormap=:magma, nan_color=:lightgray)
# sf = surface!(axT3, x3, y3, z3, color=∇T3n, colorrange=(0, 1e-4), colormap=:magma, nan_color=:lightgray)

# Colorbar(fig[2, 4], sf, height=Relative(0.6), flipaxis=true, vertical=true,
#          label="Surface temperature (ᵒC)", labelsize=20)

title = @lift @sprintf("%s days", round((times[$n]+ 1e-6) / day, digits=2))
#Label(fig[0, :], title, fontsize=24)

Label(fig[0, 1], "Latitude-Longitude Grid", tellwidth=false, fontsize=24)
Label(fig[0, 2], "Tripolar Grid", tellwidth=false, fontsize=24)
Label(fig[0, 3], "Rotated Latitude-Longitude Grid", tellwidth=false, fontsize=24)

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

save("baroclinic_wave.png", fig)

display(fig)
