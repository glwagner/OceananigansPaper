using Oceananigans
using GLMakie

filename = "flow_around_cylinder_les.jld2"
ut = FieldTimeSeries(filename, "u")
vt = FieldTimeSeries(filename, "v")

fig = Figure(size=(1000, 600))
ax = Axis(fig[2, 1], aspect=4)
slider = Slider(fig[3, 1], range=1:length(ut), startvalue=1)
n = slider.value

u = XFaceField(ut.grid)
v = YFaceField(vt.grid)

parent(u) .= parent(ut[1])
parent(v) .= parent(vt[1])

ζ = Field(∂x(v) - ∂y(u))

ζn = @lift begin
    parent(u) .= parent(ut[$n])
    parent(v) .= parent(vt[$n])
    compute!(ζ)
end

heatmap!(ax, ζn, nan_color=:gray, colormap=:balance, colorrange=(-10, 10))

t = ut.times
title = @lift string("Vorticity at t = ", t[$n])
Label(fig[1, 1], title, tellheight=true, tellwidth=false)

# record(fig, "one_flow_past_cylinder.mp4", 1:length(t), framerate=24) do nn
#     n[] = nn
# end

display(fig)

