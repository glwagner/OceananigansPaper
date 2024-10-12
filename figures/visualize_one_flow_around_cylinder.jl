using Oceananigans
using GLMakie

fig = Figure(size=(1110, 570))
ax = Axis(fig[1, 1], aspect=3, xlabel="x", ylabel="y", xaxisposition=:top)

#filename = "flow_around_cylinder_les_ReInf_Ny128_fields.jld2"
#filename = "flow_around_cylinder_les_ReInf_Ny256_fields.jld2"
#filename = "flow_around_cylinder_les_ReInf_Ny512_fields.jld2"
#filename = "flow_around_cylinder_les_ReInf_Ny1024_fields.jld2"
filename = "continued_flow_around_cylinder_les_ReInf_Ny2048_fields.jld2"

ut = FieldTimeSeries(filename, "u")
vt = FieldTimeSeries(filename, "v")

Nt = length(ut)
n = Observable(1)
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
hidedecorations!(ax)

display(fig)

Nt = length(ut)
record(fig, "flow_around_cylinder.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

