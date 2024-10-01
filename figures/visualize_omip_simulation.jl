using Oceananigans
using GLMakie

filename = "prototype_omip_simulationsurface.jld2"

ut = FieldTimeSeries(filename, "u")
et = FieldTimeSeries(filename, "e")

Nt = length(ut)
u = ut[Nt]
e = et[Nt]

fig = Figure()
axu = Axis(fig[1, 1])
axe = Axis(fig[1, 2])

heatmap!(axu, u)
heatmap!(axe, e)

display(fig)

