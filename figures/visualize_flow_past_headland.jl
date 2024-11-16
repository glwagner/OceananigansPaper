using Oceananigans
using GLMakie

ζt = FieldTimeSeries("flow_past_headland.jld2", "ζ")
ut = FieldTimeSeries("flow_past_headland.jld2", "u")
Nt = length(ut)

fig = Figure(size=(1200, 600))
axu = Axis(fig[1, 1], aspect=2)
axζ = Axis(fig[2, 1], aspect=2)
slider = Slider(fig[3, 1], range=1:Nt, startvalue=1)
n = slider.value

un = @lift interior(ut[$n], :, :, 1)
ζn = @lift interior(ζt[$n], :, :, 1)

heatmap!(axu, un)
heatmap!(axζ, ζn)

display(fig)

