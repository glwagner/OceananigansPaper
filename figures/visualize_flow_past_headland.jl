using Oceananigans
using Oceananigans.Units
using GLMakie
using Printf

Nz = 64
ut = FieldTimeSeries("flow_past_headland_$(Nz)_xy.jld2", "u")
Tt = FieldTimeSeries("flow_past_headland_$(Nz)_xy.jld2", "T")
qt = FieldTimeSeries("flow_past_headland_$(Nz)_xyz.jld2", "q")

Nq = 5 # Time index for 3D PV plot
time = qt.times[Nq]

Nu = findmin(abs.(ut.times.-time))[2] # Find corresponding time index for 2D plots

fig = Figure(size=(1100, 800))

xu, yu, zu = nodes(ut)
un = interior(ut[Nu], :, :, 1)
axu = Axis(fig[1, 1], aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", width=500)
hmu = heatmap!(axu, xu, yu, un, nan_color=:lightgray, colormap=:balance, colorrange=(-0.3, 0.3))
Colorbar(fig[0, 1], hmu, vertical=false, label="x-velocity (m s⁻¹)", width=Relative(0.5))

xT, yT, zT = nodes(Tt)
Tn = interior(Tt[Nu], :, :, 1)
axT = Axis(fig[1, 2], aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", width=500)
hmT = heatmap!(axT, xT, yT, Tn, nan_color=:lightgray, colormap=:magma, colorrange=(9.8, 11.4))
Colorbar(fig[0, 2], hmT, vertical=false, label="Temperature (ᵒC)", width=250, height=15)

H, L, δ = 256meters, 1024meters, 512meters
Nx, Ny, Nz = size(qt.grid)

# Exclude edges
xss = 2:Nx
yss = 2:Ny
zss = 2:Nz

xq, yq, zq = nodes(qt)
qt = qt[xss, yss, zss, :]
xq = xq[xss]
yq = yq[yss]
zq = zq[zss]

# Create bathymetry map
wedge(x, y) = -H *(1 + (y + abs(x)) / δ)
bathymetry = [ wedge(x, y) > z ? 1 : 0 for x=xq, y=yq, z=zq ]

ax = Axis3(fig[2, 1:2]; aspect = (ut.grid.Lx, ut.grid.Ly, 4*ut.grid.Lz), width = 1000, height=400, azimuth = 0.6π, elevation = 0.15π,
                        perspectiveness=0.8, viewmode=:fitzoom,
                        xlabel="x [m]", ylabel="y [m]", zlabel="z [m]")

# Plot bathymetry
volume!(ax, xq, yq, zq, bathymetry, algorithm = :absorption, isorange = 50, absorption=50f0, colormap = [RGBAf(0,0,0,0), :papayawhip], colorrange=(0, 1))

# adjust colormap
colormap = to_colormap(:balance)
lc = length(colormap)÷2
middle_chunk = ceil(Int, 0.4 * lc)
colormap[lc-middle_chunk:lc+middle_chunk] .= RGBAf(0, 0, 0, 0) # Make middle chunk trasparent

qₙ = qt[:, :, :, Nq]
limit = 5e-8
vol = volume!(ax, xq, yq, zq, qₙ, algorithm = :absorption, absorption=20f0, colormap=colormap, colorrange=(-limit, +limit))
Colorbar(fig, vol, bbox=ax.scene.px_area,
         label="Ertel Potential Vorticity (1/s³)", height=15, width=250, vertical=false,
         alignmode = Outside(10), halign = -0.1, valign = 0.02)

xtxt = 0.03
ytxt = 0.97
T₂ = 12.421hours
title = @sprintf("2π t / T₂ = %1.2f", time / T₂)
text!(axu, xtxt, ytxt; text=title, space=:relative, color=:black, align=(:left, :top), fontsize=18)

resize_to_layout!(fig)
save("3d_plot_headland_$Nz.png", fig)

