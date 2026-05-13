# Pre-stage NumericalEarth datasets (ETOPO bathymetry, ECCO T/S at Jan 1 1993,
# JRA55 atmosphere) on the head node before running the global-simulation
# smoke test. Datasets land in ~/.julia/datadeps/ and are visible to compute
# nodes via the shared /shared filesystem.
# Run from project root: julia --project scripts/prestage_data.jl

using Oceananigans
using Oceananigans.Units
import NumericalEarth
using Dates
using CFTime

ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"

println("=== Pre-stage: ETOPO bathymetry ===")
arch = CPU()
Nz = 40
exponential_z_faces(; Nz, depth) = -depth * (1 .- (0:Nz) / Nz).^2
z = exponential_z_faces(; Nz, depth=5000)
grid = LatitudeLongitudeGrid(arch; z,
                             size=(1440, 560, Nz),
                             halo=(7, 7, 7),
                             longitude=(0, 360),
                             latitude=(-70, 70))

bathymetry = NumericalEarth.regrid_bathymetry(grid)
println("  bathymetry: ", summary(bathymetry))

println("=== Pre-stage: ECCO T/S at 1993-01-01 ===")
dates = DateTimeProlepticGregorian(1993, 1, 1)
T_meta = NumericalEarth.ECCOMetadatum(:temperature; date=dates)
S_meta = NumericalEarth.ECCOMetadatum(:salinity; date=dates)
println("  ECCO T metadata: ", T_meta)
println("  ECCO S metadata: ", S_meta)

println("=== Pre-stage: JRA55 atmosphere ===")
atmosphere = NumericalEarth.JRA55PrescribedAtmosphere(arch)
println("  atmosphere: ", typeof(atmosphere))

println("=== Done ===")
