using Oceananigans
using CairoMakie, SixelTerm

u = FieldTimeSeries("near_global_surface_fields.jld2", "u"; backend=OnDisk())
v = FieldTimeSeries("near_global_surface_fields.jld2", "v"; backend=OnDisk())
T = FieldTimeSeries("near_global_surface_fields.jld2", "T"; backend=OnDisk())
S = FieldTimeSeries("near_global_surface_fields.jld2", "S"; backend=OnDisk())

iter = Observable(1)

times = u.times
Nt = length(times)

un = @lift(u[$iter])
vn = @lift(v[$iter])
Tn = @lift(T[$iter])
Sn = @lift(S[$iter])

uf = similar(u[1])
vf = similar(v[1])

s = Field(sqrt(uf^2 + vf^2))
ζ = Field(KernelFunctionOperation{Face, Face, Nothing}(Oceananigans.Operators.ζ₃ᶠᶠᶜ, uf.grid.underlying_grid, uf, vf))

sn = @lift begin
    set!(uf, $un)
    set!(vf, $vn)
    compute!(s)
    s
end

ζn = @lift begin
    set!(uf, $un)
    set!(vf, $vn)
    compute!(ζ)
    ζ
end

fig = Figure(size = (1200, 800))
axT = Axis(fig[2, 1], title="Surface temperature ᵒC")
axS = Axis(fig[2, 2], title="Surface salinity psu")
axs = Axis(fig[3, 1], title="Surface speed ms⁻¹")
axζ = Axis(fig[3, 2], title="Vertical vorticity s⁻¹")

hmT = heatmap!(axT, Tn, colormap=:magma, colorrange=(-1, 35))
hmS = heatmap!(axS, Sn, colormap=:haline, colorrange=(25, 35))
hms = heatmap!(axs, sn, colorrange=(0, 0.5))
hmζ = heatmap!(axζ, ζn, colorrange=(-1e-5, 1e-5), colormap=:bwr)

hidedecorations!(axT)
hidedecorations!(axS)
hidedecorations!(axs)
hidedecorations!(axζ)

Colorbar(fig[1, 1], hmT, flipaxis = false)
Colorbar(fig[1, 2], hmS, flipaxis = false)
Colorbar(fig[4, 1], hms, flipaxis = false)
Colorbar(fig[4, 2], hmζ, flipaxis = false)

iter[] = 120
CairoMakie.save("eddying_near_global.png", fig)

CairoMakie.record(fig, "eddying_near_global.mp4", 1:Nt, framerate=8) do i
    @info "Printing $i of $Nt"
    iter[] = i
end

