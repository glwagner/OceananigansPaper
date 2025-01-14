using Oceananigans
using GLMakie

#filename = "double_gyre_Nx176_Nz40_xy1.jld2"
filename = "double_gyre_Nx352_Nz40_xy1.jld2"
#filename = "double_gyre_Nx176_Nz20_xy1.jld2"
#filename = "double_gyre_Nx132_Nz30_xy1.jld2"

ut = FieldTimeSeries(filename, "u", backend=OnDisk())
vt = FieldTimeSeries(filename, "v", backend=OnDisk())
bt = FieldTimeSeries(filename, "b", backend=OnDisk())

Nt = length(ut.times)
n = Observable(Nt)

grid = ut.grid
u′ = XFaceField(grid)
v′ = YFaceField(grid)
ζ = Field(∂x(v′) - ∂y(u′))
sop = @at (Center, Center, Center) sqrt(u′^2 + v′^2)
s = Field(sop)

ζn = @lift begin
    parent(u′) .= parent(ut[$n])
    parent(v′) .= parent(vt[$n])
    compute!(ζ)
    view(ζ, :, :, grid.Nz)
end

sn = @lift begin
    parent(u′) .= parent(ut[$n])
    parent(v′) .= parent(vt[$n])
    compute!(s)
    #view(s, :, :, grid.Nz)
    interior(s, :, :, 1)
end

un = @lift interior(ut[$n], :, :, 1)
vn = @lift interior(vt[$n], :, :, 1)
bn = @lift interior(bt[$n], :, :, 1)

fig = Figure(size=(1700, 400))
axu = Axis(fig[1, 1], xlabel="Longitude (deg)", ylabel="Latitude (deg)", aspect=1.1)
axv = Axis(fig[1, 2], xlabel="Longitude (deg)", ylabel="Latitude (deg)", aspect=1.1)
axs = Axis(fig[1, 3], xlabel="Longitude (deg)", ylabel="Latitude (deg)", aspect=1.1)
axb = Axis(fig[1, 4], xlabel="Longitude (deg)", ylabel="Latitude (deg)", aspect=1.1)

bmin = minimum(bt)
bmax = maximum(bt)
Δb = bmax - bmin
ϵ = Δb/8
blims = (bmin + ϵ, bmax - ϵ)

x, y, z = nodes(bt)

heatmap!(axu, x, y, un, colormap=:balance, colorrange=(-0.5, 0.5))
heatmap!(axv, x, y, vn, colormap=:balance, colorrange=(-0.5, 0.5))
heatmap!(axs, x, y, sn, colormap=:solar, colorrange=(0, 1.0))
heatmap!(axb, x, y, bn, colormap=:magma, colorrange=blims) #(-0.31, -0.29))

display(fig)

record(fig, "double_gyre.mp4", 1:Nt, framerate=24) do nn
    @info "Plotting frame $nn of $Nt..."
    n[] = nn
end

