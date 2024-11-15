using Oceananigans
using GLMakie

xyfilename = "wind_pulse_512_xy.jld2"
yzfilename = "wind_pulse_512_yz.jld2"
xzfilename = "wind_pulse_512_xz.jld2"

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

fig = Figure(size = (900, 700))

ax = Axis3(fig[2:6, 1], aspect=(1, 1, 1/2),
           xlabel = "x (m)", ylabel = "y (m)", zlabel = "z (m)",
           xlabeloffset = 100, ylabeloffset = 100, zlabeloffset = 100,
           elevation = 0.6, azimuth = 4.1, protrusions = 40, perspectiveness = 0.7)

# slider = Slider(fig[2, 1], range=1:Nt, startvalue=Nt)
# n = slider.value
n = Observable(1)

wxy = @lift interior(wxyt[$n], :, :, 1)
wyz = @lift interior(wyzt[$n], 1, :, :)
wxz = @lift interior(wxzt[$n], :, 1, :)

wlim⁻ = Ref(1e-2)
ϵ = 0.1
wlims = @lift begin
    if $n % 10 == 0
        wmax = maximum(abs, wxyt[$n])
        wlim★ = wmax * 3/4
        wlimⁿ = ϵ * wlim⁻[] + (1 - ϵ) * wlim★
        wlim⁻[] = wlimⁿ 
    else
        wlimⁿ = wlim⁻[]
    end

    (-wlimⁿ, wlimⁿ)
end

kwargs = (colorrange=wlims, colormap=:balance)
surface!(ax, xyz, yyz, zyz; color = wyz, kwargs...)
surface!(ax, xxz, yxz, zxz; color = wxz, kwargs...)
sf = surface!(ax, xxy, yxy, zxy; color = wxy, kwargs...)

t = wxyt.times
label = @lift string("Vertical velocity (m s⁻¹) after ", prettytime(t[$n]))
Colorbar(fig[1, 1], sf; label, vertical=false, width=Relative(0.7))
rowgap!(fig.layout, Relative(-0.2))

display(fig)

record(fig, "pulse_of_wind.mp4", 1:Nt, framerate=12) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end

