using Oceananigans
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using OrthogonalSphericalShellGrids
using ClimaOcean
using GLMakie

filename = "prototype_omip_simulation_surface.jld2"

ut = FieldTimeSeries(filename, "u")
et = FieldTimeSeries(filename, "e")

Nt = length(ut)
u = interior(ut[Nt], :, :, 1)
e = interior(et[Nt], :, :, 1)

mask_immersed_field!(ut[Nt], NaN)
mask_immersed_field!(et[Nt], NaN)

fig = Figure(size=(900, 1000))
axu = Axis(fig[1, 1])
axe = Axis(fig[2, 1])

heatmap!(axu, u)
heatmap!(axe, e, colormap=:solar, colorrange=(0, 1e-3))

display(fig)

