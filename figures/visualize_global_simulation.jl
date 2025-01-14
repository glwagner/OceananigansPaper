using Oceananigans
using GLMakie

filename = "simple_global_simulation.jld2"

st = FieldTimeSeries(filename, "s")
et = FieldTimeSeries(filename, "e")
Nt = length(st)

set_theme!(Theme(fontsize=24))

fig = Figure(size=(1200, 1200))

axs = Axis(fig[1, 1])
axe = Axis(fig[2, 1])
n = Observable(1)

sn = @lift interior(st[$n], :, :, 1)
en = @lift interior(et[$n], :, :, 1)

hm = heatmap!(axs, sn, colorrange=(0, 0.7), colormap=:magma)
Colorbar(fig[1, 2], hm, label="Surface speed (m s⁻¹)")

hm = heatmap!(axe, en, colorrange=(0, 1e-3), colormap=:solar)
Colorbar(fig[2, 2], hm, label="Surface turbulent kinetic energy (m² s⁻²)")

t = st.times
title = @lift string(prettytime(t[$n]), " after Jan 1 1993")
Label(fig[0, 1:2], title, tellwidth=false)

display(fig)

record(fig, "simple_global_simulation.mp4", 1:Nt, framerate=8) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

