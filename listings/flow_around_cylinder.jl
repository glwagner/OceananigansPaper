using Oceananigans

grid = RectilinearGrid(size=(128, 512), x=(-3, 10), y=(-3, 3),
                       halo = (4, 4),
                       topology = (Periodic, Bounded, Flat))

cylinder(x, y) = (x^2 + y^2) ≤ 1
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(cylinder))

#=
@inline u_drag(x, y, z, t, u, v, w, C) = - C * u * sqrt(u^2 + v^2 + w^2)
@inline v_drag(x, y, z, t, u, v, w, C) = - C * v * sqrt(u^2 + v^2 + w^2)
@inline w_drag(x, y, z, t, u, v, w, C) = - C * w * sqrt(u^2 + v^2 + w^2)

drag_coefficient = 1e-2
u_drag_bc = FluxBoundaryCondition(u_drag, field_dependencies=(:u, :v, :w), parameters=drag_coefficient)
v_drag_bc = FluxBoundaryCondition(v_drag, field_dependencies=(:u, :v, :w), parameters=drag_coefficient)
w_drag_bc = FluxBoundaryCondition(w_drag, field_dependencies=(:u, :v, :w), parameters=drag_coefficient)

u_bcs = FieldBoundaryConditions(immersed=u_drag_bc)
v_bcs = FieldBoundaryConditions(immersed=v_drag_bc)
w_bcs = FieldBoundaryConditions(immersed=w_drag_bc)
boundary_conditions = (u=u_bcs, v=v_bcs, w=w_bcs)
=#

no_slip = ValueBoundaryCondition(0)
velocity_bcs = FieldBoundaryConditions(immersed=no_slip)
boundary_conditions = (u=velocity_bcs, v=velocity_bcs)

# DNS config
Re = 1000
advection = Centered(order=2)
closure = ScalarDiffusivity(ν=1/Re)

# LES config
advection = WENO(order=9)
closure = nothing

model = NonhydrostaticModel(; grid, closure, advection, boundary_conditions)
set!(model, u=1)

simulation = Simulation(model, Δt=1e-3, stop_time=10)

run!(simulation)

using GLMakie
u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u))
compute!(ζ)

heatmap(ζ, nan_color=:gray)
display(current_figure())

