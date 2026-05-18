# Compares three vertical mixing parameterizations against a 3D
# nonhydrostatic LES "reference" in the same single-column, wind-driven,
# surface-cooled boundary layer:
#
#   1. CATKEVerticalDiffusivity            (built-in)
#   2. TKEDissipationVerticalDiffusivity   (built-in, k-ε)
#   3. PacanowskiPhilanderVerticalDiffusivity (defined inline below — the
#      Oceananigans developer-docs reference implementation, see
#      docs/src/developer_docs/turbulence_closures.md on the Oceananigans
#      repo)
#   4. Implicit LES: 64³ NonhydrostaticModel with WENO(9) advection and no
#      explicit subgrid closure; horizontal averages serve as the truth
#      the parameterizations are trying to approximate.
#
# Produces vertical_mixing_parameterizations.png — two panels, b(z) and u(z).

using Oceananigans
using Oceananigans.Units
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures:
    CATKEVerticalDiffusivity,
    TKEDissipationVerticalDiffusivity,
    AbstractScalarDiffusivity,
    VerticalFormulation,
    VerticallyImplicitTimeDiscretization,
    buoyancy_tracers,
    buoyancy_force
using Oceananigans.BuoyancyFormulations: ∂z_b
using Oceananigans.Operators: ℑxᶜᵃᵃ, ℑyᵃᶜᵃ, ∂zᶠᶜᶠ, ∂zᶜᶠᶠ
using Oceananigans.StokesDrifts: UniformStokesDrift
using Oceananigans.Utils: launch!, prettysummary, @kernel, @index
using Statistics: mean
using Printf
using Serialization
using CUDA
using CairoMakie

include(joinpath(@__DIR__, "_smoke_prelude.jl"))

# ============================================================
# Pacanowski–Philander vertical diffusivity (Richardson-number-based).
# Reference: Pacanowski & Philander (1981), JPO. Inline implementation
# from the Oceananigans developer docs.
# ============================================================

struct PacanowskiPhilanderVerticalDiffusivity{TD, FT} <: AbstractScalarDiffusivity{TD, VerticalFormulation, 1}
    ν₀     :: FT
    ν₁     :: FT
    κ₀     :: FT
    κ_conv :: FT   # convective-adjustment diffusivity, used when N² < 0
    ν_conv :: FT   # convective-adjustment viscosity,   used when N² < 0
    c      :: FT
    n      :: FT
    maximum_diffusivity :: FT
    maximum_viscosity   :: FT
end

function PacanowskiPhilanderVerticalDiffusivity(time_discretization = VerticallyImplicitTimeDiscretization(),
                                                FT = Float64;
                                                ν₀     = 1e-3, ν₁     = 1e-1, κ₀ = 1e-4,
                                                κ_conv = 1.0,  ν_conv = 1.0,
                                                c  = 5.0,  n  = 2.0,
                                                maximum_diffusivity = Inf,
                                                maximum_viscosity   = Inf)
    TD = typeof(time_discretization)
    return PacanowskiPhilanderVerticalDiffusivity{TD, FT}(
        convert(FT, ν₀),     convert(FT, ν₁),     convert(FT, κ₀),
        convert(FT, κ_conv), convert(FT, ν_conv),
        convert(FT, c),      convert(FT, n),
        convert(FT, maximum_diffusivity), convert(FT, maximum_viscosity))
end

const PPVD = PacanowskiPhilanderVerticalDiffusivity

@inline TurbulenceClosures.viscosity_location(::PPVD)   = (Center(), Center(), Face())
@inline TurbulenceClosures.diffusivity_location(::PPVD) = (Center(), Center(), Face())
@inline TurbulenceClosures.viscosity(::PPVD, diffusivities)        = diffusivities.νz
@inline TurbulenceClosures.diffusivity(::PPVD, diffusivities, id)  = diffusivities.κz

function TurbulenceClosures.build_closure_fields(grid, clock, tracer_names, bcs, closure::PPVD)
    κz = Field{Center, Center, Face}(grid)
    νz = Field{Center, Center, Face}(grid)
    return (; κz, νz)
end

@inline ϕ²(i, j, k, grid, ϕ, args...) = ϕ(i, j, k, grid, args...)^2

@inline function shear_squaredᶜᶜᶠ(i, j, k, grid, velocities)
    ∂z_u² = ℑxᶜᵃᵃ(i, j, k, grid, ϕ², ∂zᶠᶜᶠ, velocities.u)
    ∂z_v² = ℑyᵃᶜᵃ(i, j, k, grid, ϕ², ∂zᶜᶠᶠ, velocities.v)
    return ∂z_u² + ∂z_v²
end

@inline function Riᶜᶜᶠ(i, j, k, grid, velocities, buoyancy, tracers)
    S² = shear_squaredᶜᶜᶠ(i, j, k, grid, velocities)
    N² = ∂z_b(i, j, k, grid, buoyancy, tracers)
    # 1e-10 s⁻² is the typical "minimum shear squared" floor in OGCM Ri-based
    # parameterizations — well above eps(Float64) so denom stays well-conditioned
    # for the Float64 arithmetic in ν₁/denom^(n+1).
    S²_min = convert(eltype(grid), 1e-10)
    return max(zero(grid), N²) / max(S², S²_min)
end

function TurbulenceClosures.compute_closure_fields!(diffusivities, closure::PPVD, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    tracers    = buoyancy_tracers(model)
    buoyancy   = buoyancy_force(model)
    velocities = model.velocities
    launch!(arch, grid, parameters,
            compute_pp_diffusivities!, diffusivities, grid, closure, velocities, tracers, buoyancy)
    return nothing
end

@kernel function compute_pp_diffusivities!(diffusivities, grid, closure, velocities, tracers, buoyancy)
    i, j, k = @index(Global, NTuple)
    N² = ∂z_b(i, j, k, grid, buoyancy, tracers)
    S² = shear_squaredᶜᶜᶠ(i, j, k, grid, velocities)
    S²_min = convert(eltype(grid), 1e-10)
    Ri = max(zero(grid), N²) / max(S², S²_min)
    denom = 1 + closure.c * Ri
    # Stable-stratification Pacanowski–Philander branch.
    νz_stable = closure.ν₀ + closure.ν₁ / denom^closure.n
    κz_stable = closure.κ₀ + closure.ν₀ / denom + closure.ν₁ / denom^(closure.n + 1)
    # Convective adjustment: when N² < 0, switch to (ν_conv, κ_conv).
    convecting = N² < 0
    νz = ifelse(convecting, closure.ν_conv, νz_stable)
    κz = ifelse(convecting, closure.κ_conv, κz_stable)
    νz = min(νz, closure.maximum_viscosity)
    κz = min(κz, closure.maximum_diffusivity)
    @inbounds diffusivities.νz[i, j, k] = νz
    @inbounds diffusivities.κz[i, j, k] = κz
end

Base.summary(closure::PPVD{TD}) where TD = string("PacanowskiPhilanderVerticalDiffusivity{$TD}")
function Base.show(io::IO, closure::PPVD)
    print(io, summary(closure), ":\n",
              "├── ν₀: ", prettysummary(closure.ν₀), '\n',
              "├── ν₁: ", prettysummary(closure.ν₁), '\n',
              "├── κ₀: ", prettysummary(closure.κ₀), '\n',
              "├── c:  ", prettysummary(closure.c),  '\n',
              "└── n:  ", prettysummary(closure.n))
end

# ============================================================
# Boundary-layer simulation helper.
# Same physical setup for all four cases. When `les=true`, runs a 3D
# NonhydrostaticModel with WENO(9) implicit-LES dissipation on the given
# arch. Otherwise runs the 1D HydrostaticFreeSurfaceModel with the given
# vertical-mixing closure.
# ============================================================

function boundary_layer_simulation(closure, arch=CPU(); les=false,
                                   N::Int=64,
                                   FT::DataType=Float64,
                                   N²=1e-5, Jb=1e-7, τx=-5e-4, f=1e-4,
                                   La=0.3, λ_wave=60.0,
                                   Δt=nothing,
                                   stop_time=2days)

    u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(τx))
    b_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jb))

    common_kw = (; coriolis = FPlane(; f),
                 buoyancy  = BuoyancyTracer(),
                 tracers   = :b,
                 boundary_conditions = (u = u_bcs, b = b_bcs))

    if les
        grid = RectilinearGrid(arch, FT, size=(N, N, N);
                               x=(-200, 200), y=(-200, 200), z=(-200, 0))

        # Surface waves. La_t = sqrt(u★/uˢ_surface) with u★ = sqrt(|τx|).
        # Stokes drift uˢ(z) = uˢ_surface · exp(2 k z), in +x for an
        # aligned wind/wave field. Wavelength λ_wave sets the e-folding
        # depth 1/(2k).
        if La > 0
            u★_val          = sqrt(abs(τx))
            uˢ_surface_val  = u★_val / La^2
            k_wave_val      = 2π / λ_wave
            @inline ∂z_uˢ_fn(z, t, p) = 2 * p.k * p.uˢ_surf * exp(2 * p.k * z)
            stokes_drift = UniformStokesDrift(grid; ∂z_uˢ = ∂z_uˢ_fn,
                                              parameters = (k = k_wave_val,
                                                            uˢ_surf = uˢ_surface_val))
        else
            stokes_drift = nothing
        end

        model = NonhydrostaticModel(grid; closure,
                                    advection = WENO(order=9),
                                    stokes_drift,
                                    common_kw...)
        set!(model, b = (x, y, z) -> N² * z + 1e-6 * rand())
    else
        grid = RectilinearGrid(arch, FT; size=N, z=(-200, 0),
                               topology=(Flat, Flat, Bounded))
        model = HydrostaticFreeSurfaceModel(grid; closure, common_kw...)
        set!(model, b = z -> N² * z)
    end

    # Default seed Δt:
    #   - LES: 60 s × (64/N) — at N=128 the horizontal Δx is 3.125 m and a
    #     1-minute first step blows up before the wizard's first fire (iter
    #     10) can clamp it down. The wizard grows it back from there.
    #   - 1D parameterized: 60 s, override via the Δt kwarg.
    Δt₀ = something(Δt, les ? 60.0 * 64 / N : 60.0)
    simulation = Simulation(model; Δt=Δt₀, stop_time)
    les && conjure_time_step_wizard!(simulation, cfl=0.7)

    # Progress callback — fires every 50 iters for LES (so we can see Δt and
    # max|𝐮| evolve as the boundary-layer turbulence develops) and every
    # 1000 iters for the parameterized cases (which are essentially free).
    # Plain floats only (no prettytime) — at iter 0 / startup, prettytime
    # of the Simulation can hit InexactError(Int64(-Inf)).
    wallclock = Ref(time_ns())
    prev_t    = Ref(simulation.model.clock.time)
    prev_it   = Ref(simulation.model.clock.iteration)
    function progress(sim)
        now_wall = time_ns()
        Δwall    = 1e-9 * (now_wall - wallclock[])
        Δsim_s   = sim.model.clock.time      - prev_t[]
        Δit      = sim.model.clock.iteration - prev_it[]
        spi      = Δit > 0    ? Δwall / Δit             : NaN
        sdph     = Δwall > 0  ? Δsim_s / Δwall * 3600/86400 : NaN
        u, v, w  = sim.model.velocities
        @info @sprintf("iter %d, t = %.3g h, Δt = %.3g s, wall = %.2f s, sec/iter = %.3f, SDPH ≈ %.2f, max|u,v,w| = (%.3f, %.3f, %.3f)",
                       sim.model.clock.iteration,
                       sim.model.clock.time / 3600, Float64(sim.Δt),
                       Δwall, spi, sdph,
                       maximum(abs, u), maximum(abs, v), maximum(abs, w))
        wallclock[] = now_wall
        prev_t[]    = sim.model.clock.time
        prev_it[]   = sim.model.clock.iteration
        return nothing
    end
    add_callback!(simulation, progress, IterationInterval(les ? 50 : 1000))

    smoke_test_simulation!(simulation)
    run!(simulation)
    return simulation
end

# ============================================================
# Run the four simulations: implicit-LES on GPU, three parameterizations
# on CPU (each 1D, basically free).
# ============================================================

const N_LES      = parse(Int, get(ENV, "VERT_MIX_N", "64"))
const FIG_SUFFIX = get(ENV, "VERT_MIX_FIG_SUFFIX", "")
const LES_CACHE  = get(ENV, "VERT_MIX_LES_CACHE", "les_cache$(FIG_SUFFIX).jls")

# LES profile cache: if the file exists, load it and skip the 3D run.
# Otherwise run the LES and write the file. Lets PP/CATKE/k-ε iterations
# avoid the ~30 min H100 LES once it's been computed once.
les_profile = if isfile(LES_CACHE)
    @info "Loading cached LES profile from $(LES_CACHE) (skipping 3D run)"
    deserialize(LES_CACHE)
else
    @info "Running implicit LES on GPU (3D, $(N_LES)³, WENO(9))..."
    les_arch = CUDA.functional() ? GPU() : CPU()
    les_sim = boundary_layer_simulation(nothing, les_arch; les=true, N=N_LES)
    p = let m = les_sim.model
        b̄ = compute!(Field(Average(m.tracers.b,    dims=(1, 2))))
        ū = compute!(Field(Average(m.velocities.u, dims=(1, 2))))
        (z_b = Array(znodes(m.tracers.b)),
         b   = Array(interior(b̄, 1, 1, :)),
         z_u = Array(znodes(m.velocities.u)),
         u   = Array(interior(ū, 1, 1, :)))
    end
    serialize(LES_CACHE, p)
    @info "Saved LES profile cache to $(LES_CACHE)"
    p
end

closures = ("CATKE"                => CATKEVerticalDiffusivity(),
            "k-ε"                  => TKEDissipationVerticalDiffusivity(),
            "Pacanowski–Philander" => PacanowskiPhilanderVerticalDiffusivity())

# 1D parameterized runs use their own (typically coarser) vertical grid and
# longer time-step than the LES — VerticallyImplicit closures don't have a
# diffusive CFL, so Δt is bounded only by Coriolis / explicit advection.
const N_1D = parse(Int,     get(ENV, "VERT_MIX_N_1D",       "32"))
const DT_1D = parse(Float64, get(ENV, "VERT_MIX_DT_1D_MIN", "5")) * 60

results = Dict{String, Any}()
for (name, closure) in closures
    @info "Running 1D $(name) (N = $(N_1D), Δt = $(DT_1D) s)..."
    results[name] = boundary_layer_simulation(closure, CPU(); N=N_1D, Δt=DT_1D)
end

# ============================================================
# Extract column profiles (horizontal-mean for LES, the single column
# for the parameterized cases).
# ============================================================

les_profile_tuple(p) = (p.z_b, p.b, p.z_u, p.u)

function param_profiles(sim)
    m   = sim.model
    bp  = vec(interior(m.tracers.b,    1, 1, :))
    up  = vec(interior(m.velocities.u, 1, 1, :))
    z_b = znodes(m.tracers.b)
    z_u = znodes(m.velocities.u)
    return z_b, bp, z_u, up
end

# Pull the (Face-z) tracer diffusivity from whichever field the closure exposes.
function tracer_diffusivity_field(model)
    cf = hasproperty(model, :closure_fields) ? model.closure_fields :
         hasproperty(model, :diffusivity_fields) ? model.diffusivity_fields :
         nothing
    cf === nothing && return nothing
    for k in (:κc, :κᵇ, :κu, :κz, :κ)
        hasproperty(cf, k) && return getproperty(cf, k)
    end
    return nothing
end

# ============================================================
# Plot — 3 panels: b(z), u(z), and κ(z) (log scale).
# ============================================================

set_theme!(Theme(linewidth=2.5, fontsize=14))

linestyles = Dict("CATKE" => :solid, "k-ε" => :dash, "Pacanowski–Philander" => :dashdot)
colors     = Dict("CATKE" => :steelblue, "k-ε" => :darkorange, "Pacanowski–Philander" => :seagreen)

fig = Figure(size=(1050, 350))
axb = Axis(fig[1, 1], xlabel="Buoyancy (m s⁻²)",       ylabel="z (m)")
axu = Axis(fig[1, 2], xlabel="Zonal velocity (m s⁻¹)")
axκ = Axis(fig[1, 3], xlabel="Tracer diffusivity (m² s⁻¹)",
                       xscale=log10, yaxisposition=:right, ylabel="z (m)")

z_b_les, b_les, z_u_les, u_les = les_profile_tuple(les_profile)
lines!(axb, b_les, z_b_les, label="LES", color=:black, linewidth=3)
lines!(axu, u_les, z_u_les,                color=:black, linewidth=3)
# LES has no scalar tracer diffusivity (WENO provides implicit dissipation),
# so we don't draw a black line in panel (c).

for (name, _) in closures
    sim = results[name]
    z_b, bp, z_u, up = param_profiles(sim)
    lines!(axb, bp, z_b, label=name, color=colors[name], linestyle=linestyles[name])
    lines!(axu, up, z_u,             color=colors[name], linestyle=linestyles[name])

    κ_field = tracer_diffusivity_field(sim.model)
    if κ_field !== nothing
        z_f = znodes(κ_field)
        # Drop the bottom face (impenetrable boundary, κ=0 plotted as a spike
        # on the log axis) and the top face (surface boundary, similar issue).
        # Clip the rest to avoid log10(0).
        κp  = max.(vec(interior(κ_field, 1, 1, :)), 1e-7)
        lines!(axκ, κp[2:end-1], z_f[2:end-1], color=colors[name], linestyle=linestyles[name])
    end
end

text!(axb, 0, 1; text="(a)", align=(:left, :top), space=:relative, offset=(8, -6))
text!(axu, 0, 1; text="(b)", align=(:left, :top), space=:relative, offset=(8, -6))
text!(axκ, 0, 1; text="(c)", align=(:left, :top), space=:relative, offset=(8, -6))

linkyaxes!(axb, axu, axκ)
hideydecorations!(axu, grid=false)
hidespines!(axb, :t, :r)
hidespines!(axu, :t, :r, :l)
hidespines!(axκ, :t, :l)

Legend(fig[2, 1:3], axb, orientation=:horizontal, framevisible=false)
rowsize!(fig.layout, 2, Auto(0.12))
rowgap!(fig.layout, 8)

const OUT_PNG = "vertical_mixing_parameterizations$(FIG_SUFFIX).png"
save(OUT_PNG, fig; px_per_unit=4)
@info "Saved $(OUT_PNG)"
