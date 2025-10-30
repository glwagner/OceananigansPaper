using JLD2
using Oceananigans

Nx = 2160 # longitudinal direction
Ny = 1080 # meridional direction
Nz = 1

grid = TripolarGrid(size = (Nx, Ny, Nz),
                    z = (0, 1),
                    halo = (7, 7, 7))

function geographic2cartesian(λ, φ, r=1)
    if ndims(λ) == 1
        Nλ = length(λ)
        Nφ = length(φ)
        λ = repeat(reshape(λ, Nλ, 1), 1, Nφ)
        φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)
    end

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end


λ = λnodes(grid, Center(), Center(), Center())
φ = φnodes(grid, Center(), Center(), Center())
x, y, z = geographic2cartesian(λ, φ)

filename = "stripped_checkpoint_iteration910000.jld2"

file = jldopen(filename)

@show keys(file)

uo = file["uo"]
vo = file["vo"]

ui = file["ui"]
vi = file["vi"]
ℵi = file["ℵi"]

close(file)

H = 7
Tx, Ty, Tz = size(uo)
Nx = Tx - H
Ny = Ty - H
Nz = Tz - H

uo = view(uo, 1:Nx+1, 1:Ny,   Nz)
vo = view(vo, 1:Nx,   1:Ny+1, Nz)

uo² = uo.^2
vo² = vo.^2

uo² = uo²[1:end-1, :] .+ uo²[2:end, :]
vo² = vo²[:, 1:end-1] .+ vo²[:, 2:end]
so = @. sqrt(uo² + vo²)

Tx, Ty, Tz = size(ui)
Nx = Tx - H
Ny = Ty - H
Nz = Tz - H

ui = view(ui, 1:Nx+1, 1:Ny,   1)
vi = view(vi, 1:Nx,   1:Ny+1, 1)
ℵi = view(ℵi, 1:Nx,   1:Ny, 1)

ui² = ui.^2
vi² = vi.^2

ui² = ui²[1:end-1, :] .+ ui²[2:end, :]
vi² = vi²[:, 1:end-1] .+ vi²[:, 2:end]
si = @. sqrt(ui² + vi²)

max_so = maximum(so)
max_si = maximum(si)

lim_so = max_so / 6
lim_si = max_si / 8

so[so .== 0] .= NaN
si[ℵi .< 1e-2] .= NaN

#=
fig = Figure()
axs = Axis(fig[1, 1])
axa = Axis(fig[2, 1])
heatmap!(axs, so, colorrange=(0, lim_so), colormap=:magma)
heatmap!(axs, si, colorrange=(0, lim_si), colormap=:grays)
heatmap!(axa, ℵi)
=#

fig = Figure(size=(1200, 900))

kw1 = (elevation=deg2rad(40), azimuth=deg2rad(-15), aspect=:equal)
kw2 = (elevation=deg2rad(-90), azimuth=deg2rad(165), aspect=:equal)

axs1 = Axis3(fig[1, 1]; kw1...)
axs2 = Axis3(fig[1, 2]; kw2...)

sfo = surface!(axs1, x, y, z, color=so, colorrange=(0, lim_so), colormap=:magma, nan_color=:gray90)
sfi = surface!(axs1, x, y, z, color=si, colorrange=(0, lim_si), colormap=:buda10)

sfo = surface!(axs2, x, y, z, color=so, colorrange=(0, lim_so), colormap=:magma, nan_color=:gray90)
sfi = surface!(axs2, x, y, z, color=si, colorrange=(0, lim_si), colormap=:buda10)

Colorbar(fig[2, 1], sfo, flipaxis=false, vertical=false, width=Relative(0.4),
         label="Surface ocean speed (m s⁻¹)", labelsize=16)

Colorbar(fig[2, 2], sfi, flipaxis=false, vertical=false, width=Relative(0.4),
         label="Ice speed (m s⁻¹)", labelsize=16)

for ax in (axs1, axs2)
    hidedecorations!(ax)
    hidespines!(ax)
end

colgap!(fig.layout, 1, Relative(-0.1))
rowgap!(fig.layout, 1, Relative(-0.4))

display(fig)

save("global_simulation.png", fig)
