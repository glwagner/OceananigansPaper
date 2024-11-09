using Oceananigans
using GLMakie

filename = "pulse_of_wind_xy.jld2"

wt = FieldTimeSeries(filename, "w")
Nt = length(wt)

fig = Figure()
ax = Axis(fig[1, 1])
slider = Slider(fig[2, 1], range=1:Nt, startvalue=1)
n = slider.value

wn = @lift wt[$n]

heatmap!(ax, wn)

display(fig)
