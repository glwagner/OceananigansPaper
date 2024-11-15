using Oceananigans
using GLMakie

#=
xyfilename = "eady_les_1024_xy.jld2"
yzfilename = "eady_les_1024_yz.jld2"
avfilename = "eady_les_1024_av.jld2"

ζt = FieldTimeSeries(xyfilename, "ζ")
#st = FieldTimeSeries(xyfilename, "s")
∇b²t = FieldTimeSeries(xyfilename, "∇b²")

wt = FieldTimeSeries(yzfilename, "w")
ut = FieldTimeSeries(yzfilename, "u")

Bt = FieldTimeSeries(avfilename, "b")
N²t = FieldTimeSeries(avfilename, "N²")
Ut = FieldTimeSeries(avfilename, "u")
Vt = FieldTimeSeries(avfilename, "v")
=#

Nt = length(ζt)
Nz = size(ζt, 3)

fig = Figure(size=(1800, 1200))

axζ = Axis(fig[1, 1], aspect=1, xlabel="x (m)", ylabel="y (m)")
axs = Axis(fig[1, 2], aspect=1, xlabel="x (m)", ylabel="y (m)")
axB = Axis(fig[1, 3], xlabel="Mean buoyancy (m s⁻²)", ylabel="z (m)")
axN = Axis(fig[1, 4], xlabel="Mean vertical buoyancy gradient (s⁻²)", ylabel="z (m)")
axw = Axis(fig[2, 1:4], xlabel="y (m)", ylabel="z (m)")
axu = Axis(fig[3, 1:4], xlabel="y (m)", ylabel="z (m)")

#slider = Slider(fig[4, 1:4], range=1:Nt, startvalue=1)
#n = slider.value
n = Observable(1)

ζn = @lift ζt[$n]
#sn = @lift st[$n]
∇b² = Field{Center, Center, Nothing}(∇b²t.grid)
∇b = Field(sqrt(∇b²))
∇bn = @lift begin
    parent(∇b²) .= parent(∇b²t[$n])
    compute!(∇b)
    ∇b
end
wn = @lift wt[$n]
un = @lift ut[$n]
Bn = @lift Bt[$n]
N²n = @lift N²t[$n]
Un = @lift Ut[$n]
Vn = @lift Vt[$n]

hm = heatmap!(axζ, ζn, colorrange=(-1e-3, 1e-3), colormap=:balance)
Colorbar(fig[0, 1], hm, vertical=false, label="Vertical vorticity (s⁻¹)")

hm = heatmap!(axs, ∇bn, colorrange=(0, 5e-7), colormap=:magma)
Colorbar(fig[0, 2], hm, vertical=false, label="Horizontal buoyancy gradient (s⁻²)")

hm = heatmap!(axw, wn, colorrange=(-1e-2, 1e-2), colormap=:balance)
Colorbar(fig[2, 5], hm, vertical=true, label="Vertical velocity (m s⁻¹)")

hm = heatmap!(axu, un, colorrange=(-1e-1, 1e-1), colormap=:balance)
Colorbar(fig[3, 5], hm, vertical=true, label="Zonal velocity (m s⁻¹)")

lines!(axB, Bn)
lines!(axN, N²n)

xlims!(axB, -1.2e-4, 1.2e-4)
xlims!(axN, 0, 5e-6)

colsize!(fig.layout, 3, Relative(0.1))
colsize!(fig.layout, 4, Relative(0.1))
rowsize!(fig.layout, 2, Relative(0.2))
rowsize!(fig.layout, 3, Relative(0.2))

t = ζt.times
title = @lift "Eady turbulence at " * prettytime(t[$n])
Label(fig[-1, 1:5], title)

display(fig)

record(fig, "eady_les.mp4", 1:Nt, framerate=12) do nn
    @info "Drawing frame $nn of $Nt.."
    n[] = nn
end
