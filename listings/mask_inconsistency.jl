using Oceananigans
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

grid = RectilinearGrid(size=(2, 2), x = (0, 4), z=(0, 4), topology = (Periodic, Flat, Bounded))
mask = CenterField(grid)
mask[1, 1, 1] = 1
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(mask))

c = CenterField(grid)
w = ZFaceField(grid)
set!(c, 1)
set!(w, 1)

@show maximum(c) maximum(w)

mask_immersed_field!(c, NaN)
mask_immersed_field!(w, NaN)

@show maximum(c) maximum(w)



