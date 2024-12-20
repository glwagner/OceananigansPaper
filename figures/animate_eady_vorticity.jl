using Oceananigans
using GLMakie

xyfilename = "eady_les_1024_64_xyN.jld2"
xzfilename = "eady_les_1024_64_xz1.jld2"

ζt = FieldTimeSeries(xyfilename, "ζ")
ζxzt = FieldTimeSeries(xzfilename, "ζ")

Nt = length(ζt)
Nz = size(ζt, 3)

set_theme!(Theme(fontsize=18))
fig = Figure(size=(1200, 1200))

axζ = Axis(fig[1, 1]; aspect=1)
axR = Axis(fig[2, 1]; aspect=6)

rowsize!(fig.layout, 2, Relative(0.2))

slider = Slider(fig[3, 1], startvalue=1, range=1:Nt)
n = slider.value

ζn = @lift ζt[$n]
ζxzn = @lift ζxzt[$n]

ζlims = (-1e-3, 1e-3)
ζmap = :balance

hm = heatmap!(axζ, ζn, colorrange=ζlims, colormap=ζmap)
normed = [-10, -5, 0, 5, 10]
ticks = (1e-4 .* normed, string.(normed))
Colorbar(fig[0, 1], hm; ticks, vertical=false, label="Rossby number, ζ/f",
         width = Relative(0.9))

hm = heatmap!(axR, ζxzn, colorrange=ζlims, colormap=:balance)

display(fig)

