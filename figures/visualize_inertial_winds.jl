using Oceananigans
using GLMakie

xyfilename = "pulse_of_wind_xy.jld2"
yzfilename = "pulse_of_wind_yz.jld2"
xzfilename = "pulse_of_wind_xz.jld2"

wxyt = FieldTimeSeries(xyfilename, "w")
wyzt = FieldTimeSeries(yzfilename, "w")
wxzt = FieldTimeSeries(xzfilename, "w")

Nt = length(wxyt)

x, y, z = nodes(wxzt)
Nx, Ny, Nz, Nt = size(wxzt)
j = wxzt.indices[2][1]
xxz = repeat(x, 1, Nz)
yxz = y[1] * ones(Nx, Nz)
zxz = repeat(reshape(z, 1, Nz), Nx, 1)

x, y, z = nodes(wyzt)
Nx, Ny, Nz, Nt = size(wyzt)
i = wyzt.indices[1][1]
xyz = x[i] * ones(Ny, Nz)
yyz = repeat(y, 1, Nz)
zyz = repeat(reshape(z, 1, Nz), Ny, 1)

x, y, z = nodes(wxyt)
Nx, Ny, Nz, Nt = size(wxyt)
k = wxyt.indices[3][1] + 1
xxy = x
yxy = y
zxy = z[k] * ones(Nx, Ny)

fig = Figure(size = (900, 600))

ax = Axis3(fig[1, 1], aspect=(1, 1, 1/2),
           xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)",
           xlabeloffset = 100, ylabeloffset = 100, zlabeloffset = 100,
           elevation = 0.4, azimuth = 4, protrusions = 40, perspectiveness = 0.7)

# slider = Slider(fig[2, 1], range=1:Nt, startvalue=Nt)
# n = slider.value
n = Observable(1)

# t = wxyt.times
# title = @lift "Turbulence beneath a wind pulse after " * prettytime(t[$n])
# Label(fig[1, 1], title, tellheight=false, tellwidth=false)
# rowsize!(fig.layout, 1, Relative(0.2))

wxy = @lift interior(wxyt[$n], :, :, 1)
wyz = @lift interior(wyzt[$n], 1, :, :)
wxz = @lift interior(wxzt[$n], :, 1, :)

wlim = maximum(abs, wxyt) / 16
kwargs = (; colormap=:balance)
kwargs = (colorrange=(-wlim, wlim), colormap=:balance)
surface!(ax, xyz, yyz, zyz; color = wyz, kwargs...)
surface!(ax, xxz, yxz, zxz; color = wxz, kwargs...)
surface!(ax, xxy, yxy, zxy; color = wxy, kwargs...)

#=
wxy = @lift wxyt[$n]
wyz = @lift wyzt[$n]
wxz = @lift wxzt[$n]
wlim = maximum(abs, wxyt) / 4
kwargs = (colorrange=(-wlim, wlim), colormap=:balance)
surface!(ax, wyz, :west; kwargs...)
surface!(ax, wxz, :south; kwargs...)
surface!(ax, wxy, :top; kwargs...)
=#

display(fig)

record(fig, "pulse_of_wind.mp4", 1:Nt, framerate=12) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end

