using Oceananigans
using Oceananigans.Fields: location
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using GLMakie
using Printf

Nz = 64
ut = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "u")
Tt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "T")
qt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "q")

grid = ut.grid

H, L = 256meters, 1024meters
δ = L / 2
x, y, z = (-3L, 3L), (-L, L), (-H, 0)
Nz = 64

grid = RectilinearGrid(GPU(); size=(6Nz, 2Nz, Nz), halo=(6, 6, 6),
                       x, y, z, topology=(Bounded, Bounded, Bounded))

xq, yq, zq = nodes(grid, Center(), Center(), Center())

H, L = 256meters, 1024meters
δ = L / 2

# Create bathymetry map
wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
bathymetry = [ wedge(x, y) > z ? 1 : 0 for x=xq, y=yq, z=zq ]

Nt = Observable(154)

[mask_immersed_field!(ut[n], NaN) for n in 1:length(ut.times)]
[mask_immersed_field!(Tt[n], NaN) for n in 1:length(ut.times)]

ut_plt = @lift interior(ut[$Nt], 1:200, :, :)
Tt_plt = @lift interior(Tt[$Nt], 1:200, :, :)

azimuth = -1.7π
elevation = 0.1π
perspectiveness = 0.8

fig = Figure(size=(1200, 800))

# T
ax = Axis3(fig[1, 1]; 
           aspect = (2, 1, 0.125), width = 1200, height=400, azimuth, elevation,
           perspectiveness, viewmode=:fitzoom, titlegap = 0,
           xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")
xtxt = 0.03
ytxt = 0.7
T₂ = 12.421hours
title = @lift @sprintf("t / T₂ = %1.2f", Tt.times[$Nt] / T₂)

#text!(ax, xtxt, ytxt; text=title, space=:relative, color=:black, align=(:left, :top), fontsize=18)

# Plot bathymetry
volume!(ax, (-3072, 3072), (-1024, 1024), (-256, 0), bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [:papayawhip, RGBAf(0,0,0,0), :papayawhip], colorrange=(-1, 1))

# adjust colormap

u_slices = volumeslices!(xq[1:200], yq, zq, ut_plt, colormap = Reverse(:vik), colorrange = (-0.25, 0.25), bbox_color = :transparent)

u_slices[:update_yz][](200)
u_slices[:update_xy][](grid.Nz)
u_slices[:update_xz][](grid.Ny)

xlims!(ax, -3072, 1496)
ylims!(ax, -1024, 1024)
zlims!(ax, -256, 0)

Colorbar(fig, u_slices, bbox=ax.scene.viewport,
         label="x-velocity (m/s)", height = 250, width = 15, vertical = true, halign = 0.95,
         ticks = [-0.2, -0.1, 0, 0.1, 0.2])

xq, yq, zq = nodes(Tt)

# T
ax = Axis3(fig[2, 1]; 
           aspect = (2, 1, 0.125), width = 1200, height=400, azimuth, elevation,
           perspectiveness, viewmode=:fitzoom, titlegap = 0,
           xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")

# Plot bathymetry
volume!(ax, (-3072, 3072), (-1024, 1024), (-256, 0), bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [:papayawhip, RGBAf(0,0,0,0), :papayawhip], colorrange=(-1, 1))

# adjust colormap

T_slices = volumeslices!(xq[1:200], yq, zq, Tt_plt, colormap = :lajolla, colorrange = (8.75, 12), bbox_color = :transparent)

T_slices[:update_yz][](200)
T_slices[:update_xy][](grid.Nz)
T_slices[:update_xz][](grid.Ny)


xlims!(ax, -3072, 1496)
ylims!(ax, -1024, 1024)
zlims!(ax, -256, 0)

Colorbar(fig, T_slices, bbox=ax.scene.viewport,
         label="Temperature (°C)", height = 250, width = 15, vertical = true, halign = 0.95,
         ticks = [8, 9, 10, 11, 12])

# PV
ax = Axis3(fig[3, 1]; 
           aspect = (2, 1, 0.125), width = 1200, height=400, azimuth, elevation,
           perspectiveness, viewmode=:fitzoom, titlegap = 0,
           xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")

# Plot bathymetry
volume!(ax, (-3072, 3072), (-1024, 1024), (-256, 0), bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [:papayawhip, RGBAf(0,0,0,0), :papayawhip], colorrange=(-1, 1))

# adjust colormap

colormap = to_colormap(:balance)
lc = length(colormap)÷2
middle_chunk = ceil(Int, 0.4 * lc)
colormap[lc-middle_chunk:lc+middle_chunk] .= RGBAf(0, 0, 0, 0) # Make middle chunk trasparent

qₙ = @lift interior(qt[$Nt], 2:286, 2:grid.Ny-1, 2:grid.Nz-1)
limit = 5e-8

vol = volume!(ax, (-3056, 1496), (-1008, 1008), (-252, -4), qₙ, algorithm = :absorption, absorption=20f0, colormap=colormap, colorrange=(-limit, +limit))

Colorbar(fig, vol, bbox=ax.scene.viewport,
         label="Ertel Potential Vorticity (1/s³)", height=250, width=15, vertical=true,
         halign = 0.95)

xlims!(ax, -3072, 1496)
ylims!(ax, -1024, 1024)
zlims!(ax, -256, 0)

rowgap!(fig.layout, Fixed(-100))
resize_to_layout!(fig)
save("3d_plot_headland_$(Nz)_uTq.png", fig, px_per_unit = 3)


