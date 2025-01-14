using Oceananigans
using GLMakie

xyfilename = "eady_les_1024_64_xyN.jld2"
ζt = FieldTimeSeries(xyfilename, "ζ")

Nt = length(ζt)
Nz = size(ζt, 3)

set_theme!(Theme(fontsize=18))
fig = Figure(size=(1260, 1290))

aspect = 1
xlabel = "x (km)"
ylabel = "y (km)"
ticks = [-2, -1, 0, 1, 2]
#xticks = yticks = (string.(ticks), 1e3 .* ticks)
xticks = yticks = (1e3 .* ticks, string.(ticks))
axζ1 = Axis(fig[1, 1]; aspect, xlabel, ylabel, xticks, yticks)
axζ2 = Axis(fig[1, 2]; aspect, xlabel, ylabel, xticks, yticks, yaxisposition=:right)
axζ3 = Axis(fig[2, 1]; aspect, xlabel, ylabel, xticks, yticks)
axζ4 = Axis(fig[2, 2]; aspect, xlabel, ylabel, xticks, yticks, yaxisposition=:right)

hidexdecorations!(axζ1)
hidexdecorations!(axζ2)

dx = -160
dy = 460

function xyshift(f, dx=dx, dy=dy)
    fc = circshift(interior(f), (dx, dy, 0))
    interior(f) .= fc
    return f
end

n1 = 110 + 1
ζ1 = xyshift(ζt[n1])

n2 = 145
ζ2 = xyshift(ζt[n2])

n3 = 185
ζ3 = xyshift(ζt[n3])

n4 = length(ζt)
ζ4 = xyshift(ζt[n4])

#ζlims = (-5e-4, 5e-4)
ζlims = (-1e-3, 1e-3)
ζmap = :balance

hm = heatmap!(axζ1, ζ1, colorrange=ζlims, colormap=ζmap)
normed = [-10, -5, 0, 5, 10]
ticks = (1e-4 .* normed, string.(normed))
Colorbar(fig[0, 1:2], hm; ticks, vertical=false, label="Rossby number, ζ/f",
         width = Relative(0.9))

heatmap!(axζ2, ζ2, colorrange=ζlims, colormap=ζmap)
heatmap!(axζ3, ζ3, colorrange=ζlims, colormap=ζmap)
heatmap!(axζ4, ζ4, colorrange=ζlims, colormap=ζmap)

display(fig)

@show prettytime(t[n1])
@show prettytime(t[n2])
@show prettytime(t[n3])
@show prettytime(t[n4])

x = 0.02
y = 0.98
space = :relative
align = (:left, :top)
text!(axζ1, x, y; align, space, text="4.6 days")
text!(axζ2, x, y; align, space, text="6 days")
text!(axζ3, y, x; align=(:right, :bottom), space, text="7.7 days")
text!(axζ4, x, y; align, space, text="20 days")

save("eady_vorticity.png", fig, px_per_unit=4)
