using Oceananigans
using Oceananigans.Grids: halo_size
using Oceananigans.Fields: location
using CairoMakie
using JLD2

const Nx = 4320
const Ny = 1800
const Nz = 40
const nx = 4320 ÷ 8

depth = 6000meters
z_faces = exponential_z_faces(; Nz, depth)

grid = LatitudeLongitudeGrid(CPU(), Float32;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             z = z_faces,
                             latitude  = (-75, 75),
                             longitude = (0, 360))

Hx, Hy, Hz = halo_size(grid)

function read_bathymetry(prefix)
    bottom_height = zeros(Nx, Ny)

    for rank in 0:7
        irange = nx * rank + 1 : nx * (rank + 1)
        file   = jldopen(prefix * "_$(rank).jld2")
        data   = file["serialized/grid"].immersed_boundary.bottom_height[Hx+1:nx+Hx, Hy+1:Ny+Hy, 1]
        bottom_height[irange, :] .= data
        close(file)
    end

    return bottom_height
end

bottom_height = read_bathymetry("near_global_surface_fields")

grid  = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

file0 = jldopen("near_global_surface_fields_0.jld2")
iters = keys(file0["timeseries/t"])
times = Float64[file0["timeseries/t/$(iter)"] for iter in iters]
close(file0)

utmp = FieldTimeSeries{Face,   Center, Nothing}(grid, times; backend=OnDisk(), path="near_global_surface_fields.jld2", name="u")
vtmp = FieldTimeSeries{Center, Face,   Nothing}(grid, times; backend=OnDisk(), path="near_global_surface_fields.jld2", name="v")
Ttmp = FieldTimeSeries{Center, Center, Nothing}(grid, times; backend=OnDisk(), path="near_global_surface_fields.jld2", name="T")
Stmp = FieldTimeSeries{Center, Center, Nothing}(grid, times; backend=OnDisk(), path="near_global_surface_fields.jld2", name="S")
etmp = FieldTimeSeries{Center, Center, Nothing}(grid, times; backend=OnDisk(), path="near_global_surface_fields.jld2", name="e")

function set_distributed_field_time_series!(fts, prefix)
    field = Field{location(fts)...}(grid)
    Ny = size(fts, 2)
    for (idx, iter) in enumerate(iters)
        @info "doing iter $idx of $(length(iters))"
        for rank in 0:7
            irange = nx * rank + 1 : nx * (rank + 1)
            file   = jldopen(prefix * "_$(rank).jld2")
            data   = file["timeseries/$(fts.name)/$(iter)"][Hx+1:nx+Hx, Hy+1:Ny+Hy, 1]

            interior(field, irange, :, 1) .= data
            close(file)
        end

        set!(fts, field, idx)
    end
end

set_distributed_field_time_series!(utmp, "near_global_surface_fields")
set_distributed_field_time_series!(vtmp, "near_global_surface_fields")
set_distributed_field_time_series!(Ttmp, "near_global_surface_fields")
set_distributed_field_time_series!(Stmp, "near_global_surface_fields")
set_distributed_field_time_series!(etmp, "near_global_surface_fields")
