using Oceananigans

Nx, Ny, Nz = 100, 100, 100
latitude = longitude = z = (0, 1)
underlying_grid = LatitudeLongitudeGrid(size=(Nx, Ny, Nz); latitude, longitude, z)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom((λ, φ) -> 0.5))

ci = CenterField(grid)
ciw = view(ci, 1:Nx, 1:Ny, 1:Nz)

cu = CenterField(underlying_grid)
cuw = view(cu, 1:Nx, 1:Ny, 1:Nz)

for n = 1:10
    @time minimum(ci)
    @time minimum(ciw)
    @time minimum(cu)
    @time minimum(cuw)
end
