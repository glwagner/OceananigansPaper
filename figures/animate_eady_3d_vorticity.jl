using Oceananigans
using GLMakie

xyfilename = "eady_les_1024_64_xyN.jld2"
xzfilename = "eady_les_1024_64_xz1.jld2"
yzfilename = "eady_les_1024_64_yz1.jld2"

ζxyt = FieldTimeSeries(xyfilename, "ζ")
ζxzt = FieldTimeSeries(xzfilename, "ζ")
ζyzt = FieldTimeSeries(yzfilename, "ζ")

bxyt = FieldTimeSeries(xyfilename, "∇b²")
bxzt = FieldTimeSeries(xzfilename, "∇b²")
byzt = FieldTimeSeries(yzfilename, "∇b²")

set_theme!(Theme(fontsize=18))
fig = Figure(size=(2000, 1200))

x, y, z = nodes(ζxyt)

Nx, Ny, Nz, Nt = size(ζxyt)
k = ζxyt.indices[3][1]
xt = x
yt = y
zt = z[k] * ones(Nx, Ny)

Nx, Ny, Nz, Nt = size(ζxzt)
j = ζxzt.indices[2][1]
xs = repeat(x, 1, Nz)
ys = y[j] * ones(Nx, Nz)
zs = repeat(reshape(z, 1, Nz), Nx, 1)

Nx, Ny, Nz, Nt = size(ζyzt)
i = ζyzt.indices[1][1]
xw = x[i] * ones(Ny, Nz)
yw = repeat(y, 1, Nz)
zw = repeat(reshape(z, 1, Nz), Ny, 1)

kwargs = (aspect=(1, 1, 1/16),
          xlabel = "x (m)",
          ylabel = "y (m)",
          zlabel = "z (m)",
          xlabeloffset = 100,
          ylabeloffset = 100,
          zlabeloffset = 100,
          elevation = 0.25,
          azimuth = 4.2,
          protrusions = 40,
          perspectiveness = 0.7)

axz = Axis3(fig[2, 1]; kwargs...)
axb = Axis3(fig[3, 1]; kwargs...)

#slider = Slider(fig[4, 1], startvalue=1, range=1:Nt)
#n = slider.value
n = Observable(1)

ζxyn = @lift interior(ζxyt[$n], :, :, 1)
ζxzn = @lift interior(ζxzt[$n], :, 1, :)
ζyzn = @lift interior(ζyzt[$n], 1, :, :)

bxyn = @lift sqrt.(interior(bxyt[$n], :, :, 1))
bxzn = @lift sqrt.(interior(bxzt[$n], :, 1, :))
byzn = @lift sqrt.(interior(byzt[$n], 1, :, :))

ζlims = (-1e-3, 1e-3)
ζmap = :balance

blims = (0, 1.5e-6)
bmap = :magma

surface!(axz, xt, yt, zt, color=ζxyn, colorrange=ζlims, colormap=ζmap)
surface!(axz, xs, ys, zs, color=ζxzn, colorrange=ζlims, colormap=ζmap)
hm = surface!(axz, xw, yw, zw, color=ζyzn, colorrange=ζlims, colormap=ζmap)

normed = [-10, -5, 0, 5, 10]
ticks = (1e-4 .* normed, string.(normed))
Colorbar(fig[2, 2], hm; ticks, vertical=true, height = Relative(0.4),
         label="Rossby number, ζ/f")
         

surface!(axb, xt, yt, zt, color=bxyn, colorrange=blims, colormap=bmap)
surface!(axb, xs, ys, zs, color=bxzn, colorrange=blims, colormap=bmap)
hm = surface!(axb, xw, yw, zw, color=byzn, colorrange=blims, colormap=bmap)
Colorbar(fig[3, 2], hm; vertical=true, height = Relative(0.4),
         label="Horizontal buoyancy gradient |∇ₕb| (s⁻²)")

rowgap!(fig.layout, 1, Relative(-0.1))
rowgap!(fig.layout, 2, Relative(-0.2))

t = ζxyt.times
title = @lift "Eady turbulence after t = " * prettytime(t[$n])
Label(fig[1, 1], title, tellwidth=false, tellheight=true)

display(fig)

record(fig, "eady_vorticity_buoyancy.mp4", 1:Nt, framerate=6) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end
