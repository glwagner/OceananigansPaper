using MPI

MPI.Init()
@show MPI.Comm_rank(MPI.COMM_WORLD)
@show MPI.Comm_size(MPI.COMM_WORLD)

using Oceananigans
arch = Distributed(GPU())
@show arch
