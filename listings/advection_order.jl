using Oceananigans
using Oceananigans.BoundaryConditions
using Oceananigans.TimeSteppers: compute_tendencies!

function refinement_error(scheme)
    error = []
    for N in [32, 64, 128, 256, 512, 1024]
        grid = RectilinearGrid(size=N, x=(0, 1), topology=(Periodic, Flat, Flat))  
       model = NonhydrostaticModel(grid=grid, advection=scheme, tracers=:c)

        set!(model, u = 1, c = x -> sin(6π * x))
        fill_halo_regions!(model.velocities.u, model.tracers.c)

        compute_tendencies!(model, [])
        rhs_analytical = set!(CenterField(grid), x -> - 6π * cos(6π * x)) # The analytical tendency
        rhs_numerical = model.timestepper.Gⁿ.c # The computed tendency
        push!(error, norm(rhs_analytical - rhs_numerical) / N^(0.5)) # L2 error
    end

    return error
end

errors = Dict()
orders = Dict()

label(scheme, order) = scheme == :Centered ? "Centered order $(order+1)" : string(scheme) * " order $order"
style(scheme) = scheme == :Centered ? :solid : scheme == :UpwindBiased ? :dash : :dot
colors = [:blue, :red, :green, :purple, :yellow]

for order in (1, 3, 5, 7, 9)
    error = refinement_error(Centered(order=order+1))
    errors[:Centered, order] = error
    orders[:Centered, order] = @. log10(error[1:end-1] / error[2:end]) / log10(2)
    error = refinement_error(UpwindBiased(; order))
    errors[:UpwindBiased, order] = error
    orders[:UpwindBiased, order] = @. log10(error[1:end-1] / error[2:end]) / log10(2)
    error = refinement_error(WENO(; order))
    errors[:WENO, order] = error
    orders[:WENO, order] = @. log10(error[1:end-1] / error[2:end]) / log10(2)
end

N  = [32, 64, 128, 256, 512, 1024]
plots = []

fig = Figure()
ax = Axis(fig[1, 1], xscale = log10, yscale = log10, title = "L2 errors")
for (o, order) in enumerate((1, 3, 5, 7, 9)), scheme in (:Centered, :UpwindBiased, :WENO)
    push!(plots, lines!(ax, N, errors[scheme, order], label = label(scheme, order), color = colors[o], linestyle = style(scheme)))
end
axislegend(position = :lb)

ax = Axis(fig[1, 2], xscale = log10, title = "Order")
for (o, order) in enumerate((1, 3, 5, 7, 9)), scheme in (:Centered, :UpwindBiased, :WENO)
    lines!(ax, N[2:end], orders[scheme, order], label = label(scheme, order), color = colors[o], linestyle = style(scheme))
end

plots = plots |> Array{typeof(plots[1])}
Legend(fig[1, 3], plots)
