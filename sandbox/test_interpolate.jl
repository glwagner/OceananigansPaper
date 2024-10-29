using MPI
MPI.Init()
@show MPI.has_cuda()

using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Fields: interpolate!

partition = Partition(y=2)
arch = Distributed(GPU(); partition)
x = y = z = (0, 1)
grid = RectilinearGrid(arch; size=(16, 16, 16), x, y, z, topology=(Periodic, Periodic, Bounded))
c = CenterField(grid)
set!(c, (x, y, z) -> x * y^2 * z^3)
@show c

fill_halo_regions!(c)

u = XFaceField(grid)
set!(c, (x, y, z) -> x * y^2 * z^3)
interpolate!(u, c)
@show u

