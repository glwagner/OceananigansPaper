using Oceananigans
using Oceananigans.Units
using GLMakie

filename = "random_convection_xz.jld2"
wt = FieldTimeSeries(filename, "w")
bt = FieldTimeSeries(filename, "b")
dt = FieldTimeSeries(filename, "d")

fig = Figure(size=(1400, 500))

axw = Axis(fig[1, 1], xlabel="x (m)", ylabel="z (m)", aspect=1)
axb = Axis(fig[1, 2], xlabel="x (m)", ylabel="z (m)", aspect=1)
axd = Axis(fig[1, 3], xlabel="x (m)", ylabel="z (m)", aspect=1)

Nt = length(wt)
slider = Slider(fig[2, 1:3], range=1:Nt, startvalue=1)
n = slider.value

wn = @lift wt[$n]
bn = @lift bt[$n]
dn = @lift dt[$n]

wlim = 4e-2 #maximum(abs, w) / 2

heatmap!(axw, wn, colormap=:balance, nan_color=:gray, colorrange=(-wlim, wlim))
heatmap!(axb, bn, colormap=:magma, nan_color=:gray)
heatmap!(axd, dn, colormap=:balance, nan_color=:gray)

display(fig)

