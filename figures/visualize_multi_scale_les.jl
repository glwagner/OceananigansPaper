using Oceananigans
using GLMakie

#prefix = "multi_scale_les_half_size_1024"
#prefix = "multi_scale_les_small_1024"
prefix = "multi_scale_les_big_128"
xyfilename = prefix * "_xy.jld2"
yzfilename = prefix * "_averages.jld2"

ut = FieldTimeSeries(xyfilename, "u")
vt = FieldTimeSeries(xyfilename, "v")
wt = FieldTimeSeries(xyfilename, "w")
ζt = FieldTimeSeries(xyfilename, "ζ")
bt = FieldTimeSeries(xyfilename, "b")

#uyzt = FieldTimeSeries(yzfilename, "u")
#byzt = FieldTimeSeries(yzfilename, "b")
#vyzt = FieldTimeSeries(yzfilename, "v")
#w²yzt = FieldTimeSeries(yzfilename, "w²")

t = bt.times
Nt = length(t)

fig = Figure(size=(2100, 1200))
axw = Axis(fig[1, 1], aspect=2, xlabel="x (m)", ylabel="y (m)", xaxisposition=:top)
axz = Axis(fig[2, 1], aspect=2, xlabel="x (m)", ylabel="y (m)")
axb = Axis(fig[1, 2], aspect=2, xlabel="x (m)", ylabel="y (m)", xaxisposition=:top, yaxisposition=:right)
axs = Axis(fig[2, 2], aspect=2, xlabel="x (m)", ylabel="y (m)", yaxisposition=:right)

# axeyz = Axis(fig[1, 2], aspect=1)
# axvyz = Axis(fig[2, 2], aspect=1)
# axNyz = Axis(fig[3, 2], aspect=1)

# slider = Slider(fig[3, 1:2], startvalue=1, range=1:Nt)
# n = slider.value
n = Observable(Nt)

wn = @lift wt[$n]
ζn = @lift ζt[$n]
bn = @lift bt[$n]

u1 = deepcopy(ut[1])
v1 = deepcopy(vt[1])
sop = @at (Center, Center, Center) sqrt(u1^2 + v1^2)
s = Field(sop)

sn = @lift begin
    parent(u1) .= parent(ut[$n])
    parent(v1) .= parent(vt[$n])
    compute!(s)
end

bmax = maximum(bt)
bmin = minimum(bt)
zmin = maximum(abs, ζt)
smax = 0.9 * maximum(abs, ut)

Δb = bmax - bmin
blims = (bmin + Δb/3, bmax - Δb/3)
wlim = maximum(abs, wt) / 2

heatmap!(axw, wn, colormap=:balance, colorrange=(-wlim, wlim))
heatmap!(axz, ζn, colormap=:balance, colorrange=(-5e-4, 5e-4))
heatmap!(axb, bn, colormap=:thermal, colorrange=(bmin, bmax))
heatmap!(axs, sn, colormap=:magma, colorrange=(0, smax))

#=
w²lim = maximum(w²yzt) / 4

byz = deepcopy(byzt[1])
N²yz = Field(∂z(byz))

#uyzn = @lift uyzt[$n]
eyzn = @lift w²yzt[$n]
vyzn = @lift vyzt[$n]

N²yzn = @lift begin
    parent(byz) .= parent(byzt[$n])
    compute!(N²yz)
    N²i = interior(N²yz)

    unstable = N²i .< 0
    N²i[unstable] .= NaN
    N²yz
end
=#


# heatmap!(axeyz, eyzn, colormap=:magma, colorrange=(0, w²lim))
# heatmap!(axvyz, vyzn, colormap=:balance, colorrange=(-5e-3, 5e-3))
# heatmap!(axNyz, N²yzn, colormap=:ice, colorrange=(0, 1e-5), nan_color=:red)
#contour!(axwyz, byzn, levels=15)

#display(fig)

record(fig, prefix * ".mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

