# Temporary script for visualization
using Oceananigans
using Oceananigans.Units
using GLMakie
using Printf

Nz = 32
ut = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "u")
wt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "w")
qt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "q")
ζt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "ζ")
Nt = length(qt)
x, y, z = nodes(qt)


H, L, δ = 256meters, 1024meters, 512meters
grid = qt.grid.underlying_grid

xss = 2:grid.Nx-1
yss = 2:grid.Ny-1
zss = 2:grid.Nz-1
qt = qt[xss, yss, zss, :]
x=x[xss]
y=y[yss]
z=z[zss]

wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
bathymetry = [ wedge(x, y) > z ? 1 : 0 for x=x, y=y, z=z ]

fig = Figure(size=(1100, 600))

settings_axis3 = (aspect = (grid.Lx, grid.Ly, 4*grid.Lz), azimuth = 0.6π, elevation = 0.2π,
                  perspectiveness=0.8, viewmode=:fitzoom, xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")

ax = Axis3(fig[1, 1]; settings_axis3...)
volume!(ax, x, y, z, bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [:papayawhip, RGBAf(0,0,0,0), :papayawhip], colorrange=(-1, 1))


colormap = to_colormap(:balance)
PV_lim = 5e-1
middle_chunk = ceil(Int, 0.4 * 128)
colormap[128-middle_chunk:128+middle_chunk] .= RGBAf(0, 0, 0, 0)

n = Observable(1)
qₙ = @lift qt[:,:,:,$n]
vol = volume!(ax, x, y, z, qₙ, algorithm = :absorption, absorption=20f0, colormap=colormap, colorrange=(-PV_lim, +PV_lim))
Colorbar(fig, vol, bbox=ax.scene.px_area,
         label="PV / N²∞ f₀", height=25, width=Relative(0.35), vertical=false,
         alignmode = Outside(10), halign = 0.15, valign = 0.02)


xtxt = 0.03
ytxt = 0.97
T₂ = 12.421hours
times = ζt.times
title = @lift @sprintf("2π t / T₂ = %1.2f", times[$n] / T₂)
text!(ax, xtxt, ytxt; text=title, space=:relative, color=:black, align=(:left, :top), fontsize=18)

frames = 1:4:length(times)
save("3d_plot_headland.png", fig)
