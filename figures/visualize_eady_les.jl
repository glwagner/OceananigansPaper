using Oceananigans
using GLMakie

xyfilename = "eady_les_1024_xy.jld2"
yzfilename = "eady_les_1024_yz.jld2"
avfilename = "eady_les_1024_av.jld2"

# wt = FieldTimeSeries(xyfilename, "w")
# ζt = FieldTimeSeries(xyfilename, "ζ")
# ∇b²t = FieldTimeSeries(xyfilename, "∇b²")

Nt = length(wt)
Nz = size(∇b²t, 3)

fig = Figure(size=(1300, 900))

axw1 = Axis(fig[2, 1], aspect=1, xlabel="x (m)", ylabel="y (m)")
axζ1 = Axis(fig[2, 2], aspect=1, xlabel="x (m)", ylabel="y (m)")
axb1 = Axis(fig[2, 3], aspect=1, xlabel="x (m)", ylabel="y (m)", yaxisposition=:right)
axw2 = Axis(fig[3, 1], aspect=1, xlabel="x (m)", ylabel="y (m)")
axζ2 = Axis(fig[3, 2], aspect=1, xlabel="x (m)", ylabel="y (m)")
axb2 = Axis(fig[3, 3], aspect=1, xlabel="x (m)", ylabel="y (m)", yaxisposition=:right)

hidexdecorations!(axw1)
hidexdecorations!(axζ1)
hidexdecorations!(axb1)

hideydecorations!(axζ1)
hideydecorations!(axζ2)

# slider = Slider(fig[4, 1:2], range=1:Nt, startvalue=1)
# n = slider.value

∇b² = Field{Center, Center, Nothing}(∇b²t.grid)
∇b = Field(sqrt(∇b²))

n1 = 145
w1 = wt[n1]
ζ1 = ζt[n1]

parent(∇b²) .= parent(∇b²t[n1])
compute!(∇b)
∇b1 = deepcopy(∇b)

n = Observable(633)

wn = @lift wt[$n]
ζn = @lift ζt[$n]

∇bn = @lift begin
    parent(∇b²) .= parent(∇b²t[$n])
    compute!(∇b)
    ∇b
end

wlims = (-2e-3, 2e-3)
ζlims = (-1e-3, 1e-3)
blims = (0, 5e-7)

wmap = :delta
ζmap = :balance
bmap = :magma

hm = heatmap!(axw1, w1, colorrange=wlims, colormap=wmap)
Colorbar(fig[1, 1], hm, vertical=false, label="Vertical velocity, w (m s⁻¹)",
         width = Relative(0.9))

hm = heatmap!(axζ1, ζ1, colorrange=ζlims, colormap=ζmap)
normed = [-10, -5, 0, 5, 10]
ticks = (1e-4 .* normed, string.(normed))
Colorbar(fig[1, 2], hm; ticks, vertical=false, label="Rossby number, ζ/f",
         width = Relative(0.9))

hm = heatmap!(axb1, ∇b1, colorrange=blims, colormap=bmap)
normed = [0, 1, 2, 3, 4, 5]
ticks = (1e-7 .* normed, string.(normed))
Colorbar(fig[1, 3], hm; ticks, vertical=false, width = Relative(0.9),
         label="Horizontal buoyancy gradient, |∇ₕb| (10⁻⁷ × s⁻²)")

heatmap!(axw2, wn, colorrange=wlims, colormap=wmap)
heatmap!(axζ2, ζn, colorrange=ζlims, colormap=ζmap)
heatmap!(axb2, ∇bn, colorrange=blims, colormap=bmap)

t = wt.times
title = @lift "Eady turbulence at " * prettytime(t[$n])
# Label(fig[0, 1:2], title)

display(fig)

@show prettytime(t[n1])
@show prettytime(t[n[]])

save("eady_les.png", fig)
