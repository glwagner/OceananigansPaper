using Oceananigans
using GLMakie

topology = (Bounded, Bounded, Flat)
x = y = (0, 1)
fine_grid   = RectilinearGrid(size=(1024, 1024); x, y, topology)
medium_grid = RectilinearGrid(size=(16, 16); x, y, topology)
coarse_grid = RectilinearGrid(size=(4, 4); x, y, topology)

fine_grid_y   = RectilinearGrid(size=(medium_grid.Nx, fine_grid.Ny); x, y, topology)
medium_grid_y   = RectilinearGrid(size=(coarse_grid.Nx, medium_grid.Ny); x, y, topology)

c_fine = CenterField(fine_grid)
c(x, y) = exp(x) * y
set!(c_fine, c)

c_fine_y = CenterField(fine_grid_y)
c_medium_y = CenterField(medium_grid_y)
c_medium = CenterField(medium_grid)
c_coarse = CenterField(coarse_grid)

regrid!(c_fine_y, c_fine)
regrid!(c_medium, c_fine_y)
regrid!(c_medium_y, c_medium)
regrid!(c_coarse, c_medium_y)

C1 = Field(Average(c_fine, dims=1))
C2 = Field(Average(c_medium, dims=1))
C3 = Field(Average(c_coarse, dims=1))
compute!(C1)
compute!(C2)
compute!(C3)

fig = Figure(size=(900, 300))
ax1 = Axis(fig[1, 1], aspect=1)
ax2 = Axis(fig[1, 2], aspect=1)
ax3 = Axis(fig[1, 3], aspect=1)

colormap = :magma
colorrange = (minimum(c_fine), maximum(c_fine))
heatmap!(ax1, c_fine; colorrange, colormap)
heatmap!(ax2, c_medium; colorrange, colormap)
heatmap!(ax3, c_coarse; colorrange, colormap)

hidedecorations!(ax1)
hidedecorations!(ax2)
hidedecorations!(ax3)

#=
ax4 = Axis(fig[1, 4])
y = ynodes(C1)
lines!(ax4, interior(C1, 1, :, 1), y)

y = ynodes(C2)
lines!(ax4, interior(C2, 1, :, 1), y)

y = ynodes(C3)
lines!(ax4, interior(C3, 1, :, 1), y)
hidedecorations!(ax4)
=#

display(fig)

save("finite_volumes.png", fig)
