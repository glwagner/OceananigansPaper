using Oceananigans
using Oceananigans.Fields: location
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using GLMakie
using Printf

Nz = 64
ut = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "u")
Tt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "T")

grid = ut.grid

xq, yq, zq = nodes(Tt)

H, L = 256meters, 1024meters
δ = L / 2

# Create bathymetry map
wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
bathymetry = [ wedge(x, y) > z ? 1 : 0 for x=xq, y=yq, z=zq ]

Nt = Observable(53)

[mask_immersed_field!(ut[n], NaN) for n in 1:length(ut.times)]
[mask_immersed_field!(Tt[n], NaN) for n in 1:length(ut.times)]

ut_plt = @lift interior(ut[$Nt], 1:200, :, :)
Tt_plt = @lift interior(Tt[$Nt], 1:200, :, :)

azimuth = -1.7π
elevation = 0.1π
perspectiveness = 0.8

fig = Figure(size=(1200, 800))

ax = Axis3(fig[1, 1]; 
           aspect = (2, 1, 0.125), width = 1200, height=400, azimuth, elevation,
           perspectiveness, viewmode=:fitzoom, titlegap = 0,
           xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")
xtxt = 0.03
ytxt = 0.7
T₂ = 12.421hours
title = @lift @sprintf("t / T₂ = %1.2f", Tt.times[$Nt] / T₂)

text!(ax, xtxt, ytxt; text=title, space=:relative, color=:black, align=(:left, :top), fontsize=18)

# Plot bathymetry
volume!(ax, (-3072, 3072), (-1024, 1024), (-256, 0), bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [:papayawhip, RGBAf(0,0,0,0), :papayawhip], colorrange=(-1, 1))

# adjust colormap

vols = [GLMakie.volume!(ax, (-3072, 128), (-1024, 1024), (-256, 0), ut_plt, algorithm = :iso, isovalue = u₀, isorange = 0.005, alpha = 0.95, colormap = :bam, colorrange = (-0.15, 0.15), interpolate = false) for u₀ in [-0.25:0.01:0.25;]]

xlims!(ax, -3072, 1024)

Colorbar(fig, vols[1], bbox=ax.scene.viewport,
         label="x-velocity (m/s)", height = 15, width = 250, vertical = false, valign = 0.85,
         ticks = [-0.1, 0, 0.1])

xq, yq, zq = nodes(Tt)

ax = Axis3(fig[2, 1]; 
           aspect = (2, 1, 0.125), width = 1200, height=400, azimuth, elevation,
           perspectiveness, viewmode=:fitzoom, titlegap = 0,
           xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")

# Plot bathymetry
volume!(ax, (-3072, 3072), (-1024, 1024), (-256, 0), bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [:papayawhip, RGBAf(0,0,0,0), :papayawhip], colorrange=(-1, 1))

# adjust colormap

vols = [GLMakie.volume!(ax, (-3072, 128), (-1024, 1024), (-256, 0), Tt_plt, algorithm = :iso, isovalue = T₀, isorange = 0.05, alpha = 0.95, colormap = :balance, colorrange = (7.5, 12), interpolate = false) for T₀ in [7.5:0.1:12;]]

xlims!(ax, -3072, 1024)

Colorbar(fig, vols[1], bbox=ax.scene.viewport,
         label="Temperature (°C)", height = 15, width = 250, vertical = false, valign = 0.85,
         ticks = [8, 9, 10, 11, 12])

resize_to_layout!(fig)
save("3d_plot_headland_$(Nz)_uT.png", fig, px_per_unit = 3)

