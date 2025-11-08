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

fig = Figure(size=(1400, 1200))

kw = (elevation=deg2rad(50), azimuth=deg2rad(190), aspect=:equal)

axT1 = Axis3(fig[2, 1]; kw...)
axT2 = Axis3(fig[2, 2]; kw...)

ax∇T1 = Axis3(fig[3, 1]; kw...)
ax∇T2 = Axis3(fig[3, 2]; kw...)

# slider = Slider(fig[3, 1:2], range=1:Nt, startvalue=1)
# n = slider.value
n = Observable(1)

filename1 = "baroclinic_wave_lat_lon_0"
filename2 = "baroclinic_wave_lat_lon_4"

T1 = FieldTimeSeries(filename1 * ".jld2", "T"; backend = OnDisk())
ζ1 = FieldTimeSeries(filename1 * ".jld2", "ζ"; backend = OnDisk())
∇T1 = FieldTimeSeries(filename1 * ".jld2", "∇T"; backend = OnDisk())

T2 = FieldTimeSeries(filename2 * ".jld2", "T"; backend = OnDisk())
ζ2 = FieldTimeSeries(filename2 * ".jld2", "ζ"; backend = OnDisk())
∇T2 = FieldTimeSeries(filename2 * ".jld2", "∇T"; backend = OnDisk())

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

ζ1n = @lift 1e5 * interior(ζ1[$n], :, :, 1)
T1n = @lift interior(T1[$n], :, :, 1)
∇T1n = @lift interior(∇T1[$n], :, :, 1)

ζ2n = @lift 1e5 * interior(ζ2[$n], :, :, 1)
T2n = @lift interior(T2[$n], :, :, 1)
∇T2n = @lift interior(∇T2[$n], :, :, 1)

sf = surface!(axT1, x1, y1, z1, color=T1n, colorrange=(2, 30), colormap=:thermal, nan_color=:lightgray)
sf = surface!(axT2, x2, y2, z2, color=T2n, colorrange=(2, 30), colormap=:thermal, nan_color=:lightgray)

Colorbar(fig[2, 3], sf, height=Relative(0.6), flipaxis=true, vertical=true,
         label="Surface temperature (ᵒC)", labelsize=20)

sf = surface!(ax∇T1, x1, y1, z1, color=∇T1n, colorrange=(0, 1e-4), colormap=:magma, nan_color=:lightgray)
sf = surface!(ax∇T2, x2, y2, z2, color=∇T2n, colorrange=(0, 1e-4), colormap=:magma, nan_color=:lightgray)

Colorbar(fig[3, 3], sf, height=Relative(0.6), flipaxis=true, vertical=true,
         label="Surface temperature gradient (ᵒC / m)", labelsize=20)

title = @lift @sprintf("%s days", round((times[$n]+ 1e-6) / day, digits=2))
Label(fig[0, 1:2], title, fontsize=24)

Label(fig[1, 1], "Low res", tellwidth=false, fontsize=24)
Label(fig[1, 2], "High res", tellwidth=false, fontsize=24)

for ax in (axT1, axT2, ax∇T1, ax∇T2)
    hidedecorations!(ax)
    hidespines!(ax)
end

rowgap!(fig.layout, 2, Relative(-0.1))
rowgap!(fig.layout, 3, Relative(-0.2))
colgap!(fig.layout, 1, Relative(-0.4))
colgap!(fig.layout, 2, Relative(-0.2))

display(fig)

GLMakie.record(fig, "baroclinic_resolution.mp4", 1:Nt, framerate=24) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end

