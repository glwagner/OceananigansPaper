using Oceananigans
using GLMakie

fig = Figure(size=(1110, 570))
ax1 = Axis(fig[1, 1], aspect=3, xlabel="x", ylabel="y", xaxisposition=:top)
ax2 = Axis(fig[2, 1], aspect=3, xlabel="x", ylabel="y")
ax3 = Axis(fig[3, 1], aspect=3, xlabel="x", ylabel="y")
#axF = Axis(fig[1:3, 2], xlabel="Time", ylabel="Drag coefficient")

hidexdecorations!(ax2)

n = Observable(1)

filename1 = "flow_around_cylinder_dns_Re100_Ny512_fields.jld2"
filename2 = "flow_around_cylinder_dns_Re1000_Ny2048_fields.jld2"
filename3 = "flow_around_cylinder_les_ReInf_Ny2048_fields.jld2"

ut1 = FieldTimeSeries(filename1, "u")
vt1 = FieldTimeSeries(filename1, "v")
ut2 = FieldTimeSeries(filename2, "u")
vt2 = FieldTimeSeries(filename2, "v")
ut3 = FieldTimeSeries(filename3, "u")
vt3 = FieldTimeSeries(filename3, "v")

u1 = XFaceField(ut1.grid)
u2 = XFaceField(ut2.grid)
u3 = XFaceField(ut3.grid)

v1 = YFaceField(vt1.grid)
v2 = YFaceField(vt2.grid)
v3 = YFaceField(vt3.grid)

ζ1 = Field(∂x(v1) - ∂y(u1))
ζ2 = Field(∂x(v2) - ∂y(u2))
ζ3 = Field(∂x(v3) - ∂y(u3))

ζ1n = @lift begin
    parent(u1) .= parent(ut1[$n])
    parent(v1) .= parent(vt1[$n])
    compute!(ζ1)
end

ζ2n = @lift begin
    parent(u2) .= parent(ut2[$n])
    parent(v2) .= parent(vt2[$n])
    compute!(ζ2)
end

ζ3n = @lift begin
    parent(u3) .= parent(ut3[$n])
    parent(v3) .= parent(vt3[$n])
    compute!(ζ3)
end

ζlim = 4
colorrange = (-ζlim, ζlim)
hm = heatmap!(ax1, ζ1n; nan_color=:gray, colormap=:balance, colorrange)
hm = heatmap!(ax2, ζ2n; nan_color=:gray, colormap=:balance, colorrange)
hm = heatmap!(ax3, ζ3n; nan_color=:gray, colormap=:balance, colorrange)
#Colorbar(fig[1:3, 2], hm, label="Vorticity at t=200")

labels = [
    "DNS (Re = 1)",
    "DNS (Re = 10)",
    "DNS (Re = 100)",
    "LES (Re → ∞)",
]

#=
filenames = [
    "flow_around_cylinder_dns_Re1_drag.jld2",
    "flow_around_cylinder_dns_Re10_Ny256_drag.jld2",
    "flow_around_cylinder_dns_Re100_Ny512_drag.jld2",
    "flow_around_cylinder_dns_Re1000_Ny2048_drag.jld2",
    "flow_around_cylinder_les_ReInf_Ny512_drag.jld2",
]

for (label, filename) in zip(labels, filenames)
    F = FieldTimeSeries(filename, "drag_force") 
    lines!(axF, F.times, F[1, 1, 1, :]; label)
end

Legend(fig[1:3, 3], axF)
=#

ylims!(ax1, -4, 4)
ylims!(ax2, -4, 4)
ylims!(ax3, -4, 4)

xlims!(ax1, -3, 21)
xlims!(ax2, -3, 21)
xlims!(ax3, -3, 21)

#rowsize!(fig.layout, 2, Relative(0.2))
#colsize!(fig.layout, 2, Relative(0.2))

text!(ax1, 0.01, 0.95, align=(:left, :top), space=:relative, text="DNS (Re = 50)")
text!(ax2, 0.01, 0.95, align=(:left, :top), space=:relative, text="DNS (Re = 1000)")
text!(ax3, 0.01, 0.95, align=(:left, :top), space=:relative, text="LES (Re → ∞)")

# title = @lift string("Vorticity at t = ", t[$n])
# Label(fig[1, 1], title, tellheight=true, tellwidth=false)

#hidexdecorations!(ax1)

# record(fig, "flow_past_cylinder.mp4", 1:length(t), framerate=24) do nn
#     n[] = nn
# end

n[] = 101 #length(ut1)
display(fig)

