using Oceananigans
using Oceananigans.Models.NonhydrostaticModels: ConjugateGradientPoissonSolver
using Oceananigans.Models.NonhydrostaticModels: DiagonallyDominantPreconditioner
using Oceananigans.Operators: ℑxyᶠᶜᵃ, ℑxyᶜᶠᵃ
using Oceananigans.Solvers: FFTBasedPoissonSolver
using Printf
using CUDA

@inline ϕ²(i, j, k, grid, q) = @inbounds q[i, j, k]^2
@inline speedᶠᶜᶜ(i, j, k, grid, u, v) = @inbounds sqrt(u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², v))
@inline speedᶜᶠᶜ(i, j, k, grid, u, v) = @inbounds sqrt(v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², u))
@inline u_drag(i, j, k, grid, clock, f, cᴰ) = @inbounds - cᴰ * f.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, f.u, f.v)
@inline v_drag(i, j, k, grid, clock, f, cᴰ) = @inbounds - cᴰ * f.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, f.u, f.v)

config = :les
u∞ = 1
r = 1/2
arch = GPU()
stop_time = 100

for Re in [Inf]
    if config == :dns
        if Re <= 100
            Ny = 512 
            Nx = 2Ny
        elseif Re <= 1000
            Ny = 2^11
            Nx = 2Ny
        elseif Re == 10^4
            Ny = 2^12
            Nx = 2Ny
        elseif Re == 10^5
            Ny = 2^13
            Nx = 2Ny
        elseif Re == 10^6
            Ny = 3/2 * 2^13 |> Int
            Nx = 2Ny
        end
    else
        Ny = 256
        Nx = 3Ny
        Re = Inf
    end

    cylinder(x, y) = (x^2 + y^2) ≤ r^2
    prefix = "flow_around_cylinder_$(config)_Re$(Re)_Ny$(Ny)"

    ϵ = 0 # break up-down symmetry
    x = (-6, 30) # 36
    y = (-6 + ϵ, 6 + ϵ)  # 12

    kw = (; size=(Nx, Ny), x, y, halo=(6, 6), topology=(Periodic, Bounded, Flat))
    grid = RectilinearGrid(arch; kw...)
    reduced_precision_grid = RectilinearGrid(arch, Float32; kw...)

    grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(cylinder))

    if config == :dns
        advection = Centered(order=2)
        closure = ScalarDiffusivity(ν=1/Re)

        no_slip = ValueBoundaryCondition(0)
        velocity_bcs = FieldBoundaryConditions(immersed=no_slip)
        boundary_conditions = (u=velocity_bcs, v=velocity_bcs)
    elseif config == :les
        advection = WENO(order=9)
        closure = nothing

        x₁ = minimum_xspacing(grid) / 2 
        ϰ = 0.4
        ℓ = 1e-4
        @show cᴰ = (ϰ / log(x₁ / ℓ))^2
        u_drag_bc = FluxBoundaryCondition(u_drag, discrete_form=true, parameters=cᴰ)
        v_drag_bc = FluxBoundaryCondition(v_drag, discrete_form=true, parameters=cᴰ)
        u_bcs = FieldBoundaryConditions(immersed=u_drag_bc)
        v_bcs = FieldBoundaryConditions(immersed=v_drag_bc)
        boundary_conditions = (u=u_bcs, v=v_bcs)
    end

    rate = 4
    x = xnodes(grid, Face())
    @inline mask(x, y, δ=3, x₀=27) = max(zero(x), (x - x₀) / δ)
    u_sponge = Relaxation(target=1; mask, rate)
    v_sponge = Relaxation(target=0; mask, rate)
    forcing = (u=u_sponge, v=v_sponge)

    ddp = DiagonallyDominantPreconditioner()
    preconditioner = FFTBasedPoissonSolver(reduced_precision_grid)
    reltol = abstol = 1e-7
    pressure_solver = ConjugateGradientPoissonSolver(grid, maxiter=100;
                                                     reltol, abstol, preconditioner)

    #pressure_solver = ConjugateGradientPoissonSolver(grid, maxiter=100)
    #pressure_solver = ConjugateGradientPoissonSolver(grid; preconditioner=ddp)
    #pressure_solver = nothing

    if isnothing(pressure_solver)
        prefix *= "_fft"
    end

    model = NonhydrostaticModel(; grid, pressure_solver, closure,
                                advection, forcing, boundary_conditions)

    @show model

    uᵢ(x, y) = 1e-1 * randn()
    vᵢ(x, y) = 1e-1 * randn()
    set!(model, u=uᵢ, v=vᵢ)

    Δx = minimum_xspacing(grid)
    if config == :dns
        Δt = max_Δt = 0.2 * Δx^2 * Re
    else
        Δt = 0.2 * Δx
        max_Δt = Inf
    end

    simulation = Simulation(model; Δt, stop_time)
    conjure_time_step_wizard!(simulation, cfl=1.0, IterationInterval(3); max_Δt)

    u, v, w = model.velocities
    d = ∂x(u) + ∂y(v)

    # Drag computation
    μ = XFaceField(grid)
    set!(μ, mask)
    ω = u_sponge.rate
    drag_force = Field(Integral(ω * μ * (u∞ - u)))
    compute!(drag_force)

    wall_time = Ref(time_ns())

    function progress(sim)
        if pressure_solver isa ConjugateGradientPoissonSolver
            pressure_iters = iteration(pressure_solver)
        else
            pressure_iters = 0
        end

        compute!(drag_force)
        D = CUDA.@allowscalar first(drag_force)
        cᴰ = D / (u∞ * r) 
        vmax = maximum(model.velocities.v)
        dmax = maximum(abs, d)

        msg = @sprintf("Iter: %d, time: %.2f, Δt: %.4f, Poisson iters: %d",
                       iteration(sim), time(sim), sim.Δt, pressure_iters)

        elapsed = 1e-9 * (time_ns() - wall_time[])

        msg *= @sprintf(", max d: %.2e, max v: %.2e, D: %0.2f, wall time: %s",
                        dmax, vmax, cᴰ, prettytime(elapsed))

        @info msg
        wall_time[] = time_ns()

        return nothing
    end

    add_callback!(simulation, progress, IterationInterval(100))

    ζ = ∂x(v) - ∂y(u)
    outputs = (; u, v, ζ, d)

    simulation.output_writers[:jld2] = JLD2OutputWriter(model, outputs,
                                                        schedule = TimeInterval(2),
                                                        filename = prefix * "_fields.jld2",
                                                        overwrite_existing = true,
                                                        with_halos = true)

    simulation.output_writers[:drag] = JLD2OutputWriter(model, (; drag_force),
                                                        schedule = TimeInterval(0.1),
                                                        filename = prefix * "_drag.jld2",
                                                        overwrite_existing = true,
                                                        with_halos = true)
     
    run!(simulation)
end

