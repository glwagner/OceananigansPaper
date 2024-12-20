using Oceananigans
using GLMakie

filename = "random_topography_internal_tide.jld2"
wt = FieldTimeSeries(filename, "w")
Nt = length(wt)

T₂ = 12.421hours # period of M₂ tide constituent
_, idx = findmin(abs.(wt.times .- 8T₂))

x, _, z = nodes(wt)

n = Observable(idx)
wn = @lift wt[$n]

wlim = maximum(abs, wt) * 3/4

fig = Figure(size=(1300, 200))
ax = Axis(fig[1, 1], xlabel="x (km)", ylabel="z (m)")
hm = heatmap!(ax, x ./ 1e3, z,  wn, colormap=:balance, colorrange=(-wlim, wlim), nan_color=:gray)
Colorbar(fig[1, 2], hm, label="vertical velocity (m s⁻¹)")
display(fig)

save("random_internal_tide.png", fig)

record(fig, "random_internal_tide.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
