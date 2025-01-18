using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: @compute
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using GLMakie
using Printf

Nz = 64
Tt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "T")
qt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "q")
ut = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "u")

grid = Tt.grid
xc, yc, zc = nodes(Tt)
xf, yf, zf = nodes(qt)
xf, yc, zc = nodes(ut)

H, L = 256meters, 1024meters
δ = L / 2

# Create bathymetry map
wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
bathymetry = [ wedge(x, y) > z ? 1 : 0 for x=xc, y=yc, z=zc ]

Nt = 10

mask_immersed_field!(Tt[Nt], NaN)
mask_immersed_field!(ut[Nt], NaN)

xrange_q = Nz÷3:round(Int, 4.5Nz)
qt_plt = interior(qt[Nt], xrange_q, :, :)
xrange_T = Nz÷3:3Nz
Tt_plt = interior(Tt[Nt], xrange_T, :, :)
ut_plt = interior(ut[Nt], :, :, Nz)

time = qt.times[Nt]

fig = Figure(size=(1200, 800))

ax0 = Axis(fig[0, 1]; width = 1000, height = 200, xlabel = "x [m]", ylabel = "y [m]", aspect = DataAspect())
ax0.yreversed = true
hm = heatmap!(ax0, xf, yc, ut_plt, colorrange = (-0.15, +0.15), colormap = Reverse(:RdBu), nan_color = :papayawhip)
Colorbar(fig, hm; bbox=ax0.scene.px_area, label = "u [m/s]", width = 15, tellheight = true, halign = 1.05)

axis_kwargs = (aspect = (Tt.grid.Lx, Tt.grid.Ly, 2.2*Tt.grid.Lz),
               width = 1000, height = 400, azimuth = -1.7π, elevation = 0.1π,
               perspectiveness = 0.8, viewmode = :fitzoom, titlegap = 0,
               xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]")

ax1 = Axis3(fig[1, 1]; axis_kwargs...)
ax2 = Axis3(fig[2, 1]; axis_kwargs...)

# Plot bathymetry and adjust xlim
for ax in [ax1, ax2]
    volume!(ax, xc, yc, zc, bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [RGBAf(0,0,0,0), :papayawhip], colorrange=(0, 1))
    xlims!(ax, extrema(xc[xrange_q])...)
end

cbar_kwargs = (height = 15, width = 300, vertical = false, alignmode = Outside(10), halign = 0.9, valign = 0.02)

T₂ = 12.421hours
title = @sprintf("t / T₂ = %1.2f", Tt.times[Nt] / T₂)
text!(ax1, 0.03, 0.7; text=title, space=:relative, color=:black, align=(:left, :top), fontsize=18)

plt = volumeslices!(ax1, xc[xrange_T], yc, zc, Tt_plt, colormap = Reverse(:roma), colorrange = (7.5, 12), bbox_visible = false, nan_color=:transparent)
Colorbar(fig, plt, bbox=ax1.scene.px_area, label="Temperature (°C)"; cbar_kwargs...)

plt[:update_xy][](grid.Nz-1)
plt[:update_xz][](grid.Ny-1)
plt[:update_yz][](length(xrange_T))


# adjust colormap
colormap = to_colormap(:balance)
lc = length(colormap)÷2
middle_chunk = ceil(Int, 0.4 * lc)
colormap[lc-middle_chunk:lc+middle_chunk] .= RGBAf(0, 0, 0, 0) # Make middle chunk trasparent

yrange_q = 2:grid.Ny
zrange_q = 2:grid.Nz

qt_plt = interior(qt[Nt], xrange_q, yrange_q, zrange_q)
vol = volume!(ax2, xf[xrange_q], yf[yrange_q], zf[zrange_q], qt_plt, algorithm = :absorption, absorption=20f0, colormap=colormap, colorrange=(-5e-8, +5e-8))
Colorbar(fig, vol, bbox=ax2.scene.px_area, label="Ertel Potential Vorticity (1/s³)"; cbar_kwargs...)

rowgap!(fig.layout, -50)
resize_to_layout!(fig)
save("3d_plot_headland_$Nz.png", fig)
