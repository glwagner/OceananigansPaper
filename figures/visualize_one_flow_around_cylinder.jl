using Oceananigans
using GLMakie

fig = Figure(size=(1110, 570))
ax = Axis(fig[1, 1], aspect=3, xlabel="x", ylabel="y", xaxisposition=:top)

filename = "flow_around_cylinder_les_ReInf_Ny128_fields.jld2"
#filename = "flow_around_cylinder_les_ReInf_Ny256_fields.jld2"
#filename = "flow_around_cylinder_les_ReInf_Ny512_fields.jld2"

ut = FieldTimeSeries(filename, "u")
vt = FieldTimeSeries(filename, "v")

Nt = length(ut)
slider = Slider(fig[2, 1], range=1:Nt, startvalue=1)
n = slider.value
u = XFaceField(ut.grid)
v = YFaceField(vt.grid)
ζ = Field(∂x(v) - ∂y(u))

ζn = @lift begin
    parent(u) .= parent(ut[$n])
    parent(v) .= parent(vt[$n])
    compute!(ζ)
end

ζlim = 4
colorrange = (-ζlim, ζlim)
hm = heatmap!(ax, ζn; nan_color=:gray, colormap=:balance, colorrange)
display(fig)
