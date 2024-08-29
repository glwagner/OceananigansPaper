#=
using Oceananigans
using Oceananigans.Units

# List of parameters for simulating a tide sloshing over a Gaussian mountain
Nx = 1000             # Number of cells in the horizontal direction
Nz = 200              # Number of cells in the vertical direction
Lx = 200kilometers    # m, domain width
Lz = 1000             # m, total domain height / water depth
f  = 1e-4             # s⁻¹, Coriolis parameter
hm = 300              # m, height of the Gaussian mountain
Δm = 10kilometers     # m, width of the Gaussian mountain
N² = 1e-4             # s⁻², initial buoyancy gradient
ω₂ = 2π / 12.421hours # s⁻¹, lunar semi-diurnal tidal frequency
U  = 0.01             # m s⁻¹, tidal velocity
νz = κz = 1e-4        # m² s⁻¹, vertical viscosity and diffusivity
νh = κh = 1           # m² s⁻¹, horizontal viscosity and diffusivity

# Set up a grid with a Gaussian mountain bathymetry
underlying_grid = RectilinearGrid(size = (Nx, Nz), halo = (4, 4),
                                  x = (-Lx/2, Lx/2), z = (-Lz, 0),
                                  topology = (Periodic, Flat, Bounded))

bottom_height(x) = - Lz + hm * exp(-x^2 / 2Δm^2)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

A = U * (ω₂^2 - f^2) / ω₂ # tidal amplitude
@inline tidal_forcing(x, z, t, p) = p.A * sin(p.ω₂ * t)
u_forcing = Forcing(tidal_forcing, parameters=(; ω₂, A))

# A crude model for subgrid turbulence --- or a salve for poor numerics
vertical_diffusivity   = VerticalScalarDiffusivity(ν=νz, κ=κz)
horizontal_diffusivity = HorizontalScalarDiffusivity(ν=νh, κ=κh)
closure = (vertical_diffusivity, horizontal_diffusivity)
coriolis = FPlane(; f)
buoyancy = BuoyancyTracer() # means that "b" stands for buoyancy

model = HydrostaticFreeSurfaceModel(; grid, coriolis, closure, buoyancy,
                                    tracers=:b, forcing = (; u = u_forcing))

bᵢ(x, z) = N² * z
set!(model, b=bᵢ)

simulation = Simulation(model; Δt=10, stop_time=48hours)
run!(simulation)
=#

using CairoMakie

u, v, w = model.velocities
b = model.tracers.b

B = Field(Average(b, dims=1))              # horizontally-averaged buoyancy
b′ = b - B                                 # buoyancy perturbation
q = Field(∂x(v) - ∂y(u) + f * ∂z(b′) / N²) # quasi-geostrophic potential vorticity
compute!(q)

fig = Figure(size=(1200, 400))

axw = Axis(fig[1, 1], xlabel="x (km)", ylabel="z (m)", aspect=5)
axq = Axis(fig[2, 1], xlabel="x (km)", ylabel="z (m)", aspect=5)

x, y, z = nodes(w)
heatmap!(axw, x * 1e-3, z, w, colormap=:balance, colorrange=(-1e-4, 1e-4), nan_color=:lightgray)
heatmap!(axq, x * 1e-3, z, q, colormap=:balance, colorrange=(-1e-6, 1e-6), nan_color=:lightgray)

display(fig)

