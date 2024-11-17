using Oceananigans
using GLMakie

filename = "random_topography_internal_tide.jld2"
wt = FieldTimeSeries(filename, "w")
Nt = length(wt)

n = Observable(Nt)
wn = @lift wt[$n]

wlim = maximum(abs, wt) * 3/4

fig = Figure(size=(1300, 200))
ax = Axis(fig[1, 1])
heatmap!(ax, wn, colormap=:balance, colorrange=(-wlim, wlim), nan_color=:gray)
hidedecorations!(ax)
display(fig)

save("random_internal_tide", fig)

record(fig, "internal_tide.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
