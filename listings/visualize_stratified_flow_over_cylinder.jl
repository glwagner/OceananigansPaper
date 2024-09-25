using Oceananigans
using GLMakie

filename = "stratified_flow_over_cylinder_dns_Re100_Nz256_fields.jld2"
ut = FieldTimeSeries(filename, "u")
wt = FieldTimeSeries(filename, "w")
bt = FieldTimeSeries(filename, "b")

fig = Figure(size=(1000, 600))
axu = Axis(fig[1, 1], aspect=4)
axb = Axis(fig[2, 1], aspect=4)
slider = Slider(fig[3, 1], range=1:length(ut), startvalue=1)
n = slider.value

un = @lift ut[$n]
bn = @lift bt[$n]

heatmap!(axu, un, nan_color=:gray, colormap=:magma)
heatmap!(axb, bn, nan_color=:gray, colormap=:magma)

x, y, z, = nodes(bt)
contour!(axb, x, z, bn, levels=20)
#contour!(ax, bn, levels=20)

# record(fig, "one_flow_past_cylinder.mp4", 1:length(t), framerate=24) do nn
#     n[] = nn
# end

display(fig)

