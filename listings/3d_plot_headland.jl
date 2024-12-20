# Temporary script for visualization
using Oceananigans
using Oceananigans.Units
using GLMakie
using Printf

Nz = 12
qt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "q")
Nt = length(qt)
x, y, z = nodes(qt)

H, L, δ = 256meters, 1024meters, 512meters
grid = qt.grid.underlying_grid
Nz_top = grid.Nz

wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
bathymetry = [ wedge(x, y) > z ? 1 : 0 for x=x, y=y, z=z ]

fig = Figure(size=(1100, 600))

settings_axis3 = (aspect = (grid.Lx, grid.Ly, 4*grid.Lz), azimuth = 0.6π, elevation = 0.2π,
                  perspectiveness=0.8, viewmode=:fitzoom, xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")

ax = Axis3(fig[2, 1]; settings_axis3...)
volume!(ax, x, y, z, bathymetry, algorithm = :absorption, absorption=50f0, colormap = [:papayawhip, RGBAf(0,0,0,0), :papayawhip], colorrange=(-1, 1)) # turn on anti-aliasing


colormap = to_colormap(:balance)
PV_lim = 1e-4
middle_chunk = ceil(Int, 0.5 * 128)
colormap[128-middle_chunk:128+middle_chunk] .= RGBAf(0,0,0,0)

n = Observable(1)
qₙ = @lift qt[:,:,:,$n]
vol = volume!(ax, x, y, z, qₙ, algorithm = :absorption, absorption=20f0, colormap=colormap, colorrange=(-PV_lim, +PV_lim))
Colorbar(fig, vol, bbox=ax.scene.px_area,
         label="PV / N²∞ f₀", height=25, width=Relative(0.35), vertical=false,
         alignmode = Outside(10), halign = 0.15, valign = 0.02)


pause
xtxt = 0.03
ytxt = 0.97
T₂ = 12.421hours
title = @lift @sprintf("2π t / T₂ = %1.2f", qt.times[$n] / T₂)
text!(axu, xtxt, ytxt; text=title, space=:relative, color=:white, align=(:left, :top), fontsize=18)

frames = 1:4:length(qt.times)
save("3d_plot_headland.png", fig)
