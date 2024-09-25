using Oceananigans
using GLMakie
using Printf

xzfilename = "deep_convection_no_stokes_xz.jld2"
xyfilename = "deep_convection_no_stokes_xy.jld2"

wxzt = FieldTimeSeries(xzfilename, "w")
wxyt = FieldTimeSeries(xyfilename, "w")
ζxyt = FieldTimeSeries(xyfilename, "ζ")
bxzt = FieldTimeSeries(xzfilename, "b")
bxyt = FieldTimeSeries(xyfilename, "b")
ζxzt = FieldTimeSeries(xzfilename, "ζ")

wlim = maximum(abs, wxzt) / 2
ζlim = maximum(abs, ζxzt) / 2
bmin = minimum(bxzt)
bmax = maximum(bxzt)
Δb = bmax - bmin

fig = Figure(size=(1400, 900))
axzh = Axis(fig[1, 1], aspect=1) #, aspect=4, xlabel="x (m)", ylabel="z (m)")
axwh = Axis(fig[1, 2], aspect=1) #, aspect=4, xlabel="x (m)", ylabel="z (m)")
axbh = Axis(fig[1, 3], aspect=1) #, aspect=4, xlabel="x (m)", ylabel="z (m)")
axzv = Axis(fig[2, 1], aspect=1) #, aspect=4, xlabel="x (m)", ylabel="z (m)")
axwv = Axis(fig[2, 2], aspect=1) #, aspect=4, xlabel="x (m)", ylabel="z (m)")
axbv = Axis(fig[2, 3], aspect=1) #, aspect=4, xlabel="x (m)", ylabel="z (m)")

Nt = length(wxzt)
slider = Slider(fig[3, 1:3], range=1:Nt, startvalue=1)
n = slider.value
#n = Observable(Nt)

wh = @lift wxyt[$n]
bh = @lift bxyt[$n]
ζh = @lift ζxyt[$n]

wv = @lift wxzt[$n]
bv = @lift bxzt[$n]
ζv = @lift ζxzt[$n]

heatmap!(axwv, wv, colormap=:balance, colorrange=(-wlim, wlim))
heatmap!(axwh, wh, colormap=:balance, colorrange=(-wlim, wlim))

heatmap!(axbv, bv, colormap=:magma, colorrange=(bmin + Δb/2, bmax))
heatmap!(axbh, bh, colormap=:magma, colorrange=(bmin + Δb/2, bmax))

heatmap!(axzh, ζh, colormap=:balance, colorrange=(-ζlim, ζlim))
heatmap!(axzv, ζv, colormap=:balance, colorrange=(-ζlim, ζlim))

display(fig)

# record(fig, "deep_convection.mp4", 1:Nt, framerate=24) do nn
#     @info "Drawing frame $nn of $Nt..."
#     n[] = nn
# end

