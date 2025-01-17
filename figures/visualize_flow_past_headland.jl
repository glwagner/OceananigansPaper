using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: @compute
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using GLMakie
using Printf

Nz = 64
Tt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "T")
qt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "q")

grid = Tt.grid
xT, yT, zT = nodes(Tt)

H, L = 256meters, 1024meters
δ = L / 2

# Create bathymetry map
wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
bathymetry = [ wedge(x, y) > z ? 1 : 0 for x=xT, y=yT, z=zT ]

Nt = 5

mask_immersed_field!(Tt[Nt], NaN)

@compute qt_cen = Field(@at (Center, Center, Center) qt[Nt])
xrange = 1:200
qt_plt = interior(qt_cen, xrange, :, :)
Tt_plt = interior(Tt[Nt], xrange, :, :)

axis_kwargs = (aspect = (Tt.grid.Lx, Tt.grid.Ly, 2*Tt.grid.Lz),
               width = 1000, height=400, azimuth = -1.7π, elevation = 0.1π,
               perspectiveness = 0.8, viewmode = :fitzoom, titlegap = 0,
               xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]")

cbar_kwargs = (height=15, width=250, vertical=false, alignmode = Outside(10), halign = 0.8, valign = 0.02)

time = qt.times[Nt]

fig = Figure(size=(1200, 800))

ax1 = Axis3(fig[1, 1]; axis_kwargs...)
ax2 = Axis3(fig[2, 1]; axis_kwargs...)

[ xlims!(ax, -3072, 1024) for ax in [ax1, ax2] ]

# Plot bathymetry
volume!(ax1, xT, yT, zT, bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [RGBAf(0,0,0,0), :papayawhip], colorrange=(0, 1))

T₂ = 12.421hours
title = @sprintf("t / T₂ = %1.2f", Tt.times[Nt] / T₂)
text!(ax1, 0.03, 0.7; text=title, space=:relative, color=:black, align=(:left, :top), fontsize=18)

plt = volumeslices!(ax1, xT[xrange], yT, zT, Tt_plt, colormap = Reverse(:roma), colorrange = (7.5, 12), bbox_visible = false, nan_color=:transparent)
Colorbar(fig, plt, bbox=ax1.scene.px_area,
         label="Temperature (°C)"; cbar_kwargs...)

plt[:update_xy][](grid.Nz-1)
plt[:update_xz][](grid.Ny-1)
plt[:update_yz][](xrange[end])


volume!(ax2, xq, yq, zq, bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [RGBAf(0,0,0,0), :papayawhip], colorrange=(0, 1))

# adjust colormap
colormap = to_colormap(:balance)
lc = length(colormap)÷2
middle_chunk = ceil(Int, 0.4 * lc)
colormap[lc-middle_chunk:lc+middle_chunk] .= RGBAf(0, 0, 0, 0) # Make middle chunk trasparent

Nx, Ny, Nz = size(grid)
xss = 1:300
yss = 2:Ny
zss = 2:Nz

qt_plt = interior(qt[Nt], xss, yss, zss)
vol = volume!(ax2, xq[xss], yq[yss], zq[zss], qt_plt, algorithm = :absorption, absorption=20f0, colormap=colormap, colorrange=(-5e-8, +5e-8))
Colorbar(fig, vol, bbox=ax2.scene.px_area, label="Ertel Potential Vorticity (1/s³)"; cbar_kwargs...)

resize_to_layout!(fig)
save("3d_plot_headland_$Nz.png", fig)
