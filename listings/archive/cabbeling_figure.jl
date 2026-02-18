using Oceananigans
using GLMakie
using Printf

#filename = "cabbeling_3d.jld2"
filename = "cabbeling.jld2"

ρt = FieldTimeSeries(filename, "ρ")
Tt = FieldTimeSeries(filename, "T")

Nt = length(ρt)
Nx = size(ρt, 1)

i = Int(Nx / 2)
n = Observable(121) #length(Tt))
ρ = @lift interior(ρt[$n], 1:i, 1, :)
T = @lift interior(Tt[$n], i+1:Nx, 1, :)
x, y, z = nodes(ρt)

set_theme!(Theme(fontsize=18))
fig = Figure(size=(1370, 340))

ax = Axis(fig[1, 2], aspect=4, xlabel="x (m)", ylabel="z (m)")
xlims!(ax, 0, 2)
ylims!(ax, -0.5, 0)

hm = heatmap!(ax, x[1:i], z, ρ, colormap=Makie.Reverse(:grays), colorrange=(999.90, 999.98))
Colorbar(fig[1, 1], hm, label="Density (kg m⁻³)", vertical=true, flipaxis=false)

hm = heatmap!(ax, x[i+1:end], z, T, colormap=:batlow, colorrange=(1.2, 7.3))
Colorbar(fig[1, 3], hm, label="Temperature (ᵒC)", vertical=true)

#@show title = @lift @sprintf("t = %.3f seconds", ρt.times[$n])
@show title = @sprintf("t = %.3f seconds", ρt.times[n[]])
# Label(fig[0, 1], title, tellwidth=false)

display(fig)
save("cabbeling.png", fig)

# record(fig, "cabbeling_3d.mp4", 1:Nt, framerate=24) do nn
#     @info "Drawing frame $nn of $Nt..."
#     n[] = nn
# end

