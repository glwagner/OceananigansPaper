using Oceananigans
using Oceananigans.Models.NonhydrostaticModels: ConjugateGradientPoissonSolver
using Oceananigans.Models.NonhydrostaticModels: DiagonallyDominantPreconditioner
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Solvers: FFTBasedPoissonSolver
using Printf
using CUDA
using ArgParse
using CairoMakie
using Statistics

function parse_commandline()
    s = ArgParseSettings()
  
    @add_arg_table! s begin
        "--Re"
            help = "Reynolds number"
            arg_type = Float64
            default = 5600
        "--alpha"
            help = "Slope parameter of hill"
            arg_type = Float64
            default = 1
        "--solver"
            help = "Solver to use"
            arg_type = String
            default = "cg"
    end
    return parse_args(s)
end

function floor_to_base2(x::Number)
    if x <= 0
        throw(DomainError(x, "Input must be positive"))
    end
    
    # Find the nearest power of 2 by taking log2, rounding, then raising 2 to that power
    power = floor(log2(x))
    return Int(2^power)
end

function ceil_to_base2(x::Number)
    if x <= 0
        throw(DomainError(x, "Input must be positive"))
    end
    
    # Find the nearest power of 2 by taking log2, rounding, then raising 2 to that power
    power = ceil(log2(x))
    return Int(2^power)
end

args = parse_commandline()

arch = GPU()
stop_time = 200
Re = args["Re"]
solver = args["solver"]
α = args["alpha"]

H = 1
U₀ = 1

Lx = (3.8568 * α + 5.142) * H
Ly = 0.1H
Lz = 3.036H

if Re ≈ 5600
    Nz = 512
    Nx = 1536
    Ny = 32
else
    Nz = floor_to_base2(1 / Re^-0.75 * 3.036)
    Nx = floor_to_base2(1 / Re^-0.75 * (3.8568 * α + 5.142))
    Ny = floor_to_base2(1 / Re^-0.75 * 0.1)
end

Δx = Lx / Nx
max_Δt = 0.5 * Δx^2 * Re
Δt = min(Δx / U₀ * 0.1, max_Δt)

μᵤ = 1 / (Δt * 100)

# Copied verbatim from Xiao et al. (2019), but it is wrong
# function hill(x, y, α, H)
#     x̂ = x / H / α
#     # if x̂ >= 0 && x̂ <= 0.3214α
#         # ẑ = min(1, 1 + 2.42e-4 * x̂^2 - 7.588e-5 * x̂^3)
#     # elseif x̂ > 0.3214α && x̂ <= 0.5α
#         # ẑ = 0.8955 + 3.484e-2 * x̂ - 3.629e-3 * x̂^2 + 6.749e-5 * x̂^3
#     # elseif x̂ > 0.5α && x̂ <= 0.7143α
#         # ẑ = 0.9213 + 2.931e-2 * x̂ - 3.234e-3 * x̂^2 + 5.809e-5 * x̂^3
#     # elseif x̂ > 0.7143α && x̂ <= 1.071α
#         # ẑ = 1.445 - 4.927e-2 * x̂ + 6.95e-4 * x̂^2 - 7.394e-6 * x̂^3
#     # elseif x̂ > 1.071α && x̂ <= 1.429α
#         # ẑ = 0.6401 + 3.123e-2 * x̂ - 1.988e-3 * x̂^2 + 2.242e-5 * x̂^3
#     # elseif x̂ > 1.429α && x̂ <= 1.929α
#         ẑ = max(0, 2.0139 - 7.18e-2 * x̂ + 5.875e-4 * x̂^2 + 9.553e-7 * x̂^3)
#     # else
#     #     ẑ = 0
#     # end
#     return ẑ * H
# end

function hill(x, y, z)
    x̂ = x * 28 / H / α
    if x̂ < 9
        ẑ = min(28., 28 + 6.775070969851E-03*x̂^2 - 2.1245277758E-03*x̂^3)
    elseif x̂ >= 9 && x̂ < 14
        ẑ = 2.507355893131E+01 + 9.754803562315E-01*x̂ - 1.016116352781E-01*x̂^2 + 1.889794677828E-03*x̂^3
    elseif x̂ >= 14 && x̂ < 20
        ẑ = 2.579601052357E+01 + 8.206693007457E-01*x̂ - 9.055370274339E-02*x̂^2 + 1.626510569859E-03*x̂^3
    elseif x̂ >= 20 && x̂ < 30
        ẑ = 4.046435022819E+01 - 1.379581654948E+00*x̂ + 1.945884504128E-02*x̂^2 - 2.070318932190E-04*x̂^3
    elseif x̂ >= 30 && x̂ < 40
        ẑ = 1.792461334664E+01 + 8.743920332081E-01*x̂ - 5.567361123058E-02*x̂^2 + 6.277731764683E-04*x̂^3
    elseif x̂ >= 40 && x̂ < 54
        ẑ = max(0., 5.639011190988E+01 - 2.010520359035E+00*x̂ + 1.644919857549E-02*x̂^2 + 2.674976141766E-05*x̂^3)
    else
        ẑ = 0
    end

    z_boundary = ẑ / 28 * H

    return z ≤ z_boundary
end

periodic_hill(x, y, z) = hill(x, y, z) | hill(-x + Lx, y, z)

# zs = hill.(xs, α, H) .+ hill.(-xs .+ Lx, α, H)

FILE_DIR = "./Output"
prefix = "flow_past_periodic_hill_$(solver)_Re$(Re)_alpha$(α)_Nxyz_$(Nx)x$(Ny)x$(Nz)"

x = (0, Lx)
y = (0, Ly)
z = (0, Lz)

kw = (; size=(Nx, Ny, Nz), x, y, z, halo=(6, 6, 6), topology=(Periodic, Periodic, Bounded))
grid = RectilinearGrid(arch; kw...)
reduced_precision_grid = RectilinearGrid(arch, Float32; kw...)

grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(periodic_hill))

u = XFaceField(grid)
v = YFaceField(grid)
w = ZFaceField(grid)

velocities = (; u, v, w)

u_initial(x, y, z) = U₀ + 1e-8 * randn()

set!(u, u_initial)

U_current = Field(Average(u))
compute!(U_current)
Fu = Ref(zero(grid))
Fu[] = μᵤ * (U₀ - mean(u))

auxiliary_fields = (; U_current)

@inline pressure_gradient(i, j, k, grid, clock, model_fields, Fu) = Fu[]

u_forcing = Forcing(pressure_gradient, discrete_form=true, parameters=Fu)
forcing = (; u=u_forcing)

advection = Centered(order=2)
closure = ScalarDiffusivity(ν=1/Re)
timestepper = :RungeKutta3

no_slip = ValueBoundaryCondition(0)
velocity_bcs = FieldBoundaryConditions(top=no_slip, bottom=no_slip, immersed=no_slip)
boundary_conditions = (u=velocity_bcs, v=velocity_bcs)

# ddp = DiagonallyDominantPreconditioner()
preconditioner = FFTBasedPoissonSolver(reduced_precision_grid)
reltol = abstol = 1e-7

if solver == "cg"
    pressure_solver = ConjugateGradientPoissonSolver(grid, maxiter=100;
                                                    reltol, abstol, preconditioner)
else
    pressure_solver = nothing
end

# pressure_solver = ConjugateGradientPoissonSolver(grid, maxiter=100)
# pressure_solver = ConjugateGradientPoissonSolver(grid; preconditioner=ddp)
# pressure_solver = nothing

if isnothing(pressure_solver)
    prefix *= "_fft"
end

model = NonhydrostaticModel(; grid, velocities, pressure_solver, closure,
                            advection, forcing, boundary_conditions, auxiliary_fields, timestepper)

@show model
minimum_relative_step = 1e-10

simulation = Simulation(model; Δt, stop_time, minimum_relative_step)
conjure_time_step_wizard!(simulation, cfl=0.7, IterationInterval(10); max_Δt, max_change=1.05)

wall_time = Ref(time_ns())

d = Field(∂x(u) + ∂y(v) + ∂z(w))
pNHS = model.pressures.pNHS

function progress(sim)
    if pressure_solver isa ConjugateGradientPoissonSolver
        pressure_iters = iteration(pressure_solver)
    else
        pressure_iters = 0
    end

    msg = @sprintf("Iter: %d, time: %.6f, Δt: %.4f, Poisson iters: %d",
                    iteration(sim), time(sim), sim.Δt, pressure_iters)

    elapsed = 1e-9 * (time_ns() - wall_time[])

    compute!(d)

    msg *= @sprintf(", max u: %6.3e, max w: %6.3e, max d: %6.3e, max pressure: %6.3e, Fu: %6.3e, wall time: %s",
                    maximum(abs, sim.model.velocities.u),
                    maximum(abs, sim.model.velocities.w),
                    maximum(abs, d),
                    maximum(abs, sim.model.pressures.pNHS),
                    Fu[],
                    prettytime(elapsed))

    @info msg
    wall_time[] = time_ns()

    return nothing
end

add_callback!(simulation, progress, IterationInterval(1))

function compute_bulk_statistics(sim)
    compute!(sim.model.auxiliary_fields.U_current)
    # Fu[] += μᵤ * (U₀ - mean(sim.model.auxiliary_fields.U_current))
    Fu[] = μᵤ * (U₀ - mean(sim.model.auxiliary_fields.U_current))
end

add_callback!(simulation, compute_bulk_statistics, IterationInterval(1))

OUTPUT_DIR = "$(FILE_DIR)/$(prefix)"
mkpath(OUTPUT_DIR)
@info "Output directory: $OUTPUT_DIR"

ζ = ∂z(u) - ∂x(w)

ubar = Average(u, dims=2)
vbar = Average(v, dims=2)
wbar = Average(w, dims=2)
ζbar = Average(ζ, dims=2)
dbar = Average(d, dims=2)
pNHSbar = Average(pNHS, dims=2)

outputs = (; ubar, vbar, wbar, ζbar, dbar, pNHSbar)

simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                    schedule = TimeInterval(1),
                                                    filename = "$(OUTPUT_DIR)/fields.jld2",
                                                    with_halos = true)

simulation.output_writers[:averaged_jld2] = JLD2OutputWriter(model, outputs,
                                                             filename = "$(OUTPUT_DIR)/averaged_fields.jld2",
                                                             schedule = AveragedTimeInterval(20, window=20),
                                                             with_halos = true)

simulation.output_writers[:checkpointer] = Checkpointer(model,
                                                        schedule = TimeInterval(40),
                                                        prefix = "$(OUTPUT_DIR)/checkpointer")

run!(simulation)

# xF = xnodes(grid, Face())
# xC = xnodes(grid, Center())
# yC = ynodes(grid, Center())
# zC = znodes(grid, Center())

# #%%
# fig = Figure()
# ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y", aspect=DataAspect())
# CUDA.@allowscalar hm = heatmap!(ax, xF, zC, Array(interior(u, :, 1, :)))
# # hm = heatmap!(ax, xC, zC, interior(grid.immersed_boundary.mask, :, 1, :))
# Colorbar(fig[1, 2], hm, tellheight=false)
# display(fig)
#%%
ubar_data = FieldTimeSeries("$(OUTPUT_DIR)/fields.jld2", "ubar")
vbar_data = FieldTimeSeries("$(OUTPUT_DIR)/fields.jld2", "vbar")
wbar_data = FieldTimeSeries("$(OUTPUT_DIR)/fields.jld2", "wbar")

xF = xnodes(ubar_data.grid, Face())
xC = xnodes(ubar_data.grid, Center())
yC = ynodes(ubar_data.grid, Center())
yF = ynodes(ubar_data.grid, Face())
zC = znodes(ubar_data.grid, Center())
zF = znodes(ubar_data.grid, Face())
Nt = length(ubar_data.times)

ulim = (-maximum(abs, ubar_data), maximum(abs, ubar_data))
vlim = (-maximum(abs, vbar_data), maximum(abs, vbar_data))
wlim = (-maximum(abs, wbar_data), maximum(abs, wbar_data))
#%%

fig = Figure(size=(800, 1200))
n = Observable(1)
axu = Axis(fig[1, 1], xlabel = "x", ylabel = "z", aspect=DataAspect(), title = "u")
axv = Axis(fig[2, 1], xlabel = "x", ylabel = "z", aspect=DataAspect(), title = "v")
axw = Axis(fig[3, 1], xlabel = "x", ylabel = "z", aspect=DataAspect(), title = "w")

uₙ = @lift interior(ubar_data[$n], :, 1, :)
vₙ = @lift interior(vbar_data[$n], :, 1, :)
wₙ = @lift interior(wbar_data[$n], :, 1, :)

hmu = heatmap!(axu, xF, zC, uₙ, colormap = :balance, colorrange = ulim)
hmv = heatmap!(axv, xC, zC, vₙ, colormap = :balance, colorrange = vlim)
hmw = heatmap!(axw, xC, zF, wₙ, colormap = :balance, colorrange = wlim)

Colorbar(fig[1, 2], hmu)
Colorbar(fig[2, 2], hmv)
Colorbar(fig[3, 2], hmw)

title = @lift "Re = $(Re), t = $(ubar_data.times[$n])"
Label(fig[0, :], title, font=:bold, tellwidth=false)
trim!(fig.layout)

CairoMakie.record(fig, "./$(OUTPUT_DIR)/$(prefix)_velocities.mp4", 1:Nt, framerate=3, px_per_unit=2) do nn
    @info "Recording frame $nn"
    n[] = nn
end
#%%
ζbar_data = FieldTimeSeries("$(OUTPUT_DIR)/fields.jld2", "ζbar")

xF = xnodes(ζbar_data.grid, Face())
xC = xnodes(ζbar_data.grid, Center())
yC = ynodes(ζbar_data.grid, Center())
yF = ynodes(ζbar_data.grid, Face())
zC = znodes(ζbar_data.grid, Center())
zF = znodes(ζbar_data.grid, Face())
Nt = length(ζbar_data.times)

ζlim = (-maximum(abs, ζbar_data), maximum(abs, ζbar_data))
#%%
fig = Figure(size=(800, 500))
n = Observable(1)
axζ = Axis(fig[1, 1], xlabel = "x", ylabel = "z", aspect=DataAspect(), title = "ζ")

ζₙ = @lift interior(ζbar_data[$n], :, 1, :)

hmζ = heatmap!(axζ, xF, zF, ζₙ, colormap = :balance, colorrange = ζlim)

Colorbar(fig[1, 2], hmζ)

title = @lift "Re = $(Re), t = $(ubar_data.times[$n])"
Label(fig[0, :], title, font=:bold, tellwidth=false)
trim!(fig.layout)

CairoMakie.record(fig, "./$(OUTPUT_DIR)/$(prefix)_vorticity.mp4", 1:Nt, framerate=3, px_per_unit=2) do nn
    @info "Recording frame $nn"
    n[] = nn
end
#%%