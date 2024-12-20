# Temporary script for visualization
using Oceananigans
using Oceananigans.Units
using CairoMakie
using Printf

Nz = 8#64
ut = FieldTimeSeries("flow_past_headland_$(Nz)_xy.jld2", "u")
Tt = FieldTimeSeries("flow_past_headland_$(Nz)_xy.jld2", "T")
wt = FieldTimeSeries("flow_past_headland_$(Nz)_xy.jld2", "w")
ζt = FieldTimeSeries("flow_past_headland_$(Nz)_xy.jld2", "ζ")
Nt = length(ut)
x, y, z = nodes(ut)

Nz_top = ut.grid.underlying_grid.Nz

fig = Figure(size=(1100, 600))
n = Observable(1)

un = @lift interior(ut[$n], :, :, 1)
axu = Axis(fig[1, 1], aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", width=500)
hmu = heatmap!(axu, x, y, un, nan_color=:lightgray, colormap=:balance, colorrange=(-0.3, 0.3))
Colorbar(fig[0, 1], hmu, vertical=false, label="x-velocity (m s⁻¹)", width=Relative(0.5))

Tn = @lift interior(Tt[$n], :, :, 1)
axT = Axis(fig[1, 2], aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", width=500)
hmT = heatmap!(axT, x, y, Tn, nan_color=:lightgray, colormap=:magma, colorrange=(9.8, 11.4))
Colorbar(fig[0, 2], hmT, vertical=false, label="Temperature (ᵒC)", width=Relative(0.5))

wn = @lift interior(wt[$n], :, :, 1)
axw = Axis(fig[3, 1], aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", width=500)
hmw = heatmap!(axw, x, y, wn, nan_color=:lightgray, colormap=:balance, colorrange=(-0.03, 0.03))
Colorbar(fig[2, 1], hmw, vertical=false, label="z-velocity (m s⁻¹)", width=Relative(0.5))

ζn = @lift interior(ζt[$n], :, :, 1)
axζ = Axis(fig[3, 2], aspect=DataAspect(), xlabel="x (m)", ylabel="y (m)", width=500)
hmζ = heatmap!(axζ, x, y, ζn, nan_color=:lightgray, colormap=:balance, colorrange=(-0.005, 0.005))
Colorbar(fig[2, 2], hmζ, vertical=false, label="z-vorticity (s⁻¹)", width=Relative(0.5))

xtxt = 0.03
ytxt = 0.97
T₂ = 12.421hours
title = @lift @sprintf("2π t / T₂ = %1.2f", ut.times[$n] / T₂)
text!(axu, xtxt, ytxt; text=title, space=:relative, color=:white, align=(:left, :top), fontsize=18)

frames = 1:length(ut.times)
record(fig, "flow_past_headland_$Nz.mp4", frames, framerate=14) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end
