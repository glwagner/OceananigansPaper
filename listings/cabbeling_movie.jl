using Oceananigans
using GLMakie
using Printf

filename = "cabbeling_3d.jld2"
#filename = "cabbeling.jld2"

ρt = FieldTimeSeries(filename, "ρ")
Tt = FieldTimeSeries(filename, "T")

Nt = length(ρt)
Nx = size(ρt, 1)
i = Int(Nx / 2)
n = Observable(1)
ρ = @lift interior(ρt[$n], 1:i, 1, :)
T = @lift interior(Tt[$n], i+1:Nx, 1, :)
x, y, z = nodes(ρt)

fig = Figure(size=(1050, 750))

ax = Axis(fig[2, 1], aspect=2, xlabel="x (m)", ylabel="z (m)")
xlims!(ax, 0.0, 1.0)
ylims!(ax, -0.5, 0)

hm = heatmap!(ax, x[1:i], z, ρ, colormap=Makie.Reverse(:grays), colorrange=(999.90, 999.98))
Colorbar(fig[1, 1], hm, label="Density (kg m⁻³)", vertical=false)

hm = heatmap!(ax, x[i+1:end], z, T, colormap=:magma, colorrange=(1.2, 7.3))
Colorbar(fig[3, 1], hm, label="Temperature (ᵒC)", vertical=false, flipaxis=false)

title = @lift @sprintf("t = %.3f seconds", ρt.times[$n])
Label(fig[0, 1], title, tellwidth=false)

display(fig)

record(fig, "cabbeling_3d.mp4", 1:Nt, framerate=24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

