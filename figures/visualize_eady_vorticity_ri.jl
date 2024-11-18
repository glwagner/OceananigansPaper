using Oceananigans
using GLMakie

xyNfilename = "eady_les_1024_64_xyN.jld2"
xy2filename = "eady_les_1024_64_xy2.jld2"

# yzfilename = "eady_les_1024_64_yz.jld2"
# avfilename = "eady_les_1024_64_av.jld2"

Rt = FieldTimeSeries(xy2filename, "Ri", backend=OnDisk())
ζt = FieldTimeSeries(xyNfilename, "ζ", backend=OnDisk())

Nt = length(ζt)
Nz = size(ζt, 3)

fig = Figure(size=(1200, 1200))

axζ1 = Axis(fig[2, 1], aspect=1, xlabel="x (m)", ylabel="y (m)")
axR1 = Axis(fig[2, 2], aspect=1, xlabel="x (m)", ylabel="y (m)", yaxisposition=:right)

axζ2 = Axis(fig[3, 1], aspect=1, xlabel="x (m)", ylabel="y (m)")
axR2 = Axis(fig[3, 2], aspect=1, xlabel="x (m)", ylabel="y (m)", yaxisposition=:right)

hidexdecorations!(axζ1)
hidexdecorations!(axR1)

# slider = Slider(fig[4, 1:2], range=1:Nt, startvalue=1)
# n = slider.value

R = Field{Center, Center, Nothing}(Rt.grid)
R⁻¹ = Field(1 / R)

n1 = 145
ζ1 = ζt[n1]

parent(R) .= parent(Rt[n1])
compute!(R⁻¹)
R1 = deepcopy(R⁻¹)

n = Observable(481)

Rc = zeros(size(parent(Rt[1]))...)
ζc = zeros(size(parent(ζt[1]))...)

dx = -150
dy = 450

Rn = @lift begin
    circshift!(parent(R), parent(Rt[$n]), (dx, dy, 0))
    compute!(R⁻¹)
    R⁻¹
end

ζn = @lift begin
    circshift!(ζc, parent(ζt[$n]), (dx, dy, 0))
    parent(ζt[$n]) .= ζc
    ζt[$n]
end

ζlims = (-1e-3, 1e-3)
Rlims = (0.0, 5)

wmap = :delta
Rmap = :magma

hm = heatmap!(axζ1, ζ1, colorrange=ζlims, colormap=ζmap)
normed = [-10, -5, 0, 5, 10]
ticks = (1e-4 .* normed, string.(normed))
Colorbar(fig[1, 1], hm; ticks, vertical=false, label="Rossby number, ζ/f",
         width = Relative(0.9))

hm = heatmap!(axR1, R1, colorrange=Rlims, colormap=Rmap)
normed = [0, 1, 2, 3, 4, 5]
ticks = (1e-7 .* normed, string.(normed))
Colorbar(fig[1, 2], hm; vertical=false, width = Relative(0.9),
         label="Richardson number")

heatmap!(axζ2, ζn, colorrange=ζlims, colormap=ζmap)
heatmap!(axR2, Rn, colorrange=Rlims, colormap=Rmap)

display(fig)

@show prettytime(t[n1])
@show prettytime(t[n[]])

#save("eady_les.png", fig)
