using Oceananigans
using GLMakie

ℵ  = FieldTimeSeries("sea_ice_surface_fields.jld2", "ℵ")
ui = FieldTimeSeries("sea_ice_surface_fields.jld2", "ui")
vi = FieldTimeSeries("sea_ice_surface_fields.jld2", "vi")
uo = FieldTimeSeries("ocean_surface_fields.jld2", "uo")
vo = FieldTimeSeries("ocean_surface_fields.jld2", "vo")

n = Observable(1)

u1 = similar(ui[1])
v1 = similar(vi[1])
ℵ1 = similar(ℵ[1])

bat = interior(uo[400], :, :, 1) .== 0
bat = Float64.(bat)
bat[bat .== 1] .= NaN;

u2 = similar(uo[1])
v2 = similar(vo[1])

si = Field(sqrt(u1^2 + v1^2))
so = Field(sqrt(u2^2 + v2^2))

sin = @lift begin
    set!(u1, ui[$n])
    set!(v1, vi[$n])
    Oceananigans.BoundaryConditions.fill_halo_regions!((u1, v1, ℵ1, u2, v2))
    compute!(si)
    sii = interior(si, :, :, 1) .* (interior(ℵ[$n], :, :, 1) .> 1e-3)
    sii[sii .== 0.0] .= NaN
    sii
end

son = @lift begin
    set!(u2, uo[$n * 2])
    set!(v2, vo[$n * 2])
    Oceananigans.BoundaryConditions.fill_halo_regions!((u2, v2))
    compute!(so)
    interior(so, :, :, 1) .+ bat
end

fig = Figure(size = (1200, 600))
ax1 = Axis(fig[1, 1:2])
hm1 = heatmap!(ax1, son; colorrange = (0.0, 0.8), colormap = :deep, nan_color=:grey)
hm2 = heatmap!(ax1, sin; colorrange = (0.0, 0.5), colormap = Reverse(:greys))
hidespines!(ax1)
hidedecorations!(ax1)
Colorbar(fig[2, 1], hm1, vertical=false, ticks = ([0.2, 0.4, 0.6], [L"0.2", L"0.4", L"0.6"]), flipaxis=false, label=L"\text{\textbf{Ocean speed (m/s)}}")
Colorbar(fig[2, 2], hm2, vertical=false, ticks = ([0.2, 0.4], [L"0.2", L"0.4"]), flipaxis=false, label=L"\text{\textbf{Seaice speed (m/s)}}")
