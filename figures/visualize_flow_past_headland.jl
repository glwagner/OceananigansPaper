using Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using CairoMakie
using Printf
#using MathTeXEngine

ut = FieldTimeSeries("flow_past_headland_64_xy.jld2", "u")
Tt = FieldTimeSeries("flow_past_headland_64_xy.jld2", "T")
Nt = length(ut)
x, y, z = nodes(ut)

T₂ = 12.421hours

# times
tp = T₂ .* (1 .+ [1/3, 2/3, 1])#2.75, 5.25]
np = map(t -> findfirst(ut.times .>= t), tp)

# mask_immersed_field! is giving me a segfault.
# This is a workaround.
Nx, Ny, Nz, Nt = size(ut)
Nz_top = ut.grid.underlying_grid.Nz
for n in np
    @inbounds for i = 1:Nx, j=1:Ny
        if ut[i, j, Nz_top, n] == 0
            ut[i, j, Nz_top, n] = NaN
        end

        if Tt[i, j, Nz_top, n] == 0
            Tt[i, j, Nz_top, n] = NaN
        end
    end

    # We shouldn't have to do this manually either
    # mask_immersed_field!(ut[n], NaN)
    # mask_immersed_field!(Tt[n], NaN)
end

fig = Figure(size=(1100, 900))

for (m, n) in enumerate(np)
    axu = Axis(fig[m, 1], xlabel="x (m)", ylabel="y (m)", aspect = DataAspect())
    axz = Axis(fig[m, 2], xlabel="x (m)", ylabel="y (m)", yaxisposition=:right, aspect = DataAspect())

    xlims!(axu, -2000, 2000)
    xlims!(axz, -2000, 2000)

    un = interior(ut[n], :, :, 1)
    Tn = interior(Tt[n], :, :, 1)

    hmu = heatmap!(axu, x, y, un, nan_color=:lightgray, colormap=:balance, colorrange=(-0.3, 0.3))
    hmT = heatmap!(axz, x, y, Tn, nan_color=:lightgray, colormap=:magma, colorrange=(8.7, 11.4))

    if m == 1
        Colorbar(fig[0, 1], hmu, vertical=false, label="x-velocity (m s⁻¹)", width=Relative(0.9))
        Colorbar(fig[0, 2], hmT, vertical=false, label="Temperature (ᵒC)", width=Relative(0.9))
    end

    if m < 3
        hidexdecorations!(axu)
        hidexdecorations!(axz)
    end

    xtxt = 0.03
    ytxt = 0.97
    t = ut.times[n]/T₂
    label = @sprintf("t / T₂ = %.2f", t)
    if m == 1
        color = :white
    elseif m == 2
        color = :black
    elseif m == 3
        color = :black
    end
    text!(axu, xtxt, ytxt; text=label, space=:relative, color, align=(:left, :top), fontsize=18)
end

display(fig)

save("flow_past_headland.png", fig)

