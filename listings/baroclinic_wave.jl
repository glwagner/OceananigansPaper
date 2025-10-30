using Oceananigans
using Oceananigans.OrthogonalSphericalShellGrids: RotatedLatitudeLongitudeGrid
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Printf
using CUDA

# for grid_type in ("lat_lon", "rotated_lat_lon", "tripolar")

# grid_type = "lat_lon"
# grid_type = "rotated_lat_lon"
# grid_type = "tripolar"
# grid_type = "cubed_sphere"

arch = GPU()
#resolution = 1
resolution = 2
Nx = 360 ÷ resolution
Ny = 170 ÷ resolution
Nz = 10
H = 3000
halo = (7, 7, 7)
latitude = (-85, 85)
longitude = (0, 360)
z = (-H, 0)
size = (Nx, Ny, Nz)

if resolution == 2
    filename = "baroclinic_wave_" * grid_type * "_0"
else
    filename = "baroclinic_wave_" * grid_type * "_" * string(Int(1/resolution))
end

grid = if grid_type == "lat_lon"
  LatitudeLongitudeGrid(arch; size, halo, latitude, longitude, z)
elseif grid_type == "rotated_lat_lon"
    RotatedLatitudeLongitudeGrid(arch; size, halo, latitude, longitude, z, north_pole=(70, 55))
elseif grid_type == "tripolar"
    # TripolarGrid with "Gaussian islands" over the two north poles
    dφ, dλ = 4, 8
    λ₀, φ₀ = 70, 55
    h = 100
    
    isle(λ, φ) = exp(-(λ - λ₀)^2 / 2dλ^2 - (φ - φ₀)^2 / 2dφ^2)
    gaussian_isles(λ, φ) = - H + (H + h) * (isle(λ, φ) + isle(λ - 180, φ))

    cylinder(λ, φ) = ((λ - λ₀)^2 / 2dλ^2 + (φ - φ₀)^2 / 2dφ^2) < 1
    cylindrical_isles(λ, φ) = - H + (H + h) * (cylinder(λ, φ) + cylinder(λ - 180, φ))
    
    underlying_grid = TripolarGrid(arch; size, halo, z)
    #ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(gaussian_isles))
    ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(cylindrical_isles))
elseif grid_type == "cubed_sphere"
    Np = Int(Nx / 6)
    ConformalCubedSphereGrid(arch, panel_size=(Np, Np, Nz),
                             z_halo=7, horizontal_direction_halo=7, z=(-H, 0))
else
    error("Dont understand grid_type $grid_type")
end
                                             
momentum_advection = WENOVectorInvariant(order=9)
tracer_advection = WENO(order=7)
coriolis = HydrostaticSphericalCoriolis()
buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())
free_surface = SplitExplicitFreeSurface(grid, substeps=60)

model = HydrostaticFreeSurfaceModel(; grid, coriolis, free_surface,
                                      buoyancy, tracers=(:T, :S),
                                      momentum_advection, tracer_advection)

Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 45) / 8)) / 2 + rand()
Sᵢ(λ, φ, z) = 28 - 5e-3 * z + rand()
set!(model, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=10minutes, stop_time=180days)

function progress(sim)
    T = sim.model.tracers.T
    S = sim.model.tracers.S
    u, v, w = sim.model.velocities

    msg = @sprintf("%d: %s, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim),
                   maximum(abs, u),
                   maximum(abs, v),
                   maximum(abs, w))

    msg *= @sprintf(", extrema(T): (%.2f, %.2f), extrema(S): (%.2f, %.2f)",
                    minimum(T), maximum(T), minimum(S), maximum(S))

    @info msg

    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
ζ = ∂x(v) - ∂y(u)
∇T = sqrt(∂x(T)^2 + ∂y(T)^2)
#fields = (; u, v, w, T, S, ζ)
fields = (; ζ, T, ∇T)

ow = JLD2Writer(model, fields,
                filename = filename * ".jld2",
                indices = (:, :, grid.Nz),
                schedule = TimeInterval(0.5days),
                overwrite_existing = true)

simulation.output_writers[:surface] = ow

run!(simulation)

