using Oceananigans
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using GLMakie

xy_filename = "double_gyre_Nx704_Nz40_xy1.jld2"
xz_filename = "double_gyre_Nx704_Nz40_xz.jld2"

uxy = FieldTimeSeries(xy_filename, "u", backend=OnDisk())
vxy = FieldTimeSeries(xy_filename, "v", backend=OnDisk())
bxy = FieldTimeSeries(xy_filename, "b", backend=OnDisk())

uxz = FieldTimeSeries(xz_filename, "u", backend=OnDisk())
vxz = FieldTimeSeries(xz_filename, "v", backend=OnDisk())
bxz = FieldTimeSeries(xz_filename, "b", backend=OnDisk())

times = bxy.times
Nt = length(times)
grid = bxy.grid
λ, φ, z = nodes(bxy)

Nx, Ny, Nz = size(grid)

λ_xz = repeat(λ, 1, Nz)
φ_xz = φ[1] * ones(Nx, Nz)
z_xz = repeat(reshape(z, 1, Nz), Nx, 1)

λ_xy = λ
φ_xy = φ
z_xy = z[end] * ones(Nx, Ny)

bnxy = interior(bxy[end], :, :, 1)
bnxz = interior(bxz[end], :, 1, :)

unxy = uxy[end]
vnxy = vxy[end]
ζnxy = Field(∂x(vnxy) - ∂y(unxy))
compute!(ζnxy)

sop = @at (Center, Center, Center) sqrt(unxy^2 + vnxy^2)
snxy = Field(sop)
compute!(snxy)

#vnxy = interior(vxy[end], :, :, 1)
#vnxz = interior(vxz[end], :, 1, :)

#mask_immersed_field!(bxy[end], NaN)
#mask_immersed_field!(bxz[end], NaN)

fig = Figure(size = (1600, 800))

kwargs = (aspect = (1, 1, 0.2),
          xlabel = "Longitude (deg)",
          ylabel = "Latitude (deg)",
          zlabel = "z (m)",
          xlabeloffset = 100,
          ylabeloffset = 100,
          zlabeloffset = 100,
          elevation = 0.45,
          azimuth = 6.8,
          xspinesvisible = false,
          zgridvisible = false,
          protrusions = 40,
          perspectiveness = 0.7)

#axb = Axis3(fig[1, 1]; kwargs...)
#axu = Axis3(fig[1, 2]; kwargs...)
# surface!(axb, λ_xz, φ_xz, z_xz; color=bnxz, colorrange=(-1e-2, 0))
# surface!(axb, λ_xy, φ_xy, z_xy; color=bnxy, colorrange=(-1e-2, 0))
# surface!(axu, λ_xz, φ_xz, z_xz; color=vnxz)
# surface!(axu, λ_xy, φ_xy, z_xy; color=vnxy)
#surface!(axu, λ_xy, φ_xy, z_xy; color=ζnxy, colormap=:balance, colorrange=(-1e-4, 1e-4))

axs = Axis(fig[1, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)")
axu = Axis(fig[1, 2], xlabel="Longitude (deg)", ylabel="Latitude (deg)")

#mask_immersed_field!(snxy, NaN)
#mask_immersed_field!(ζnxy, NaN)

heatmap!(axs, λ_xy, φ_xy, interior(snxy, :, :, 1), colormap=:magma, colorrange=(0, 1))
heatmap!(axu, λ_xy, φ_xy, interior(ζnxy, :, :, 1), colormap=:balance, colorrange=(-2e-4, 2e-4))

display(fig)

