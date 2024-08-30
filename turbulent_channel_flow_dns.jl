using Oceananigans

grid = RectilinearGrid(GPU(),
                       size=(512, 256, 256),
                       x = (0, 4π)
                       y = (-π, π)
                       z = (-1, 1)
                       topology = (Periodic, Bounded, Bounded))

no_slip = ValueBoundaryCondition(0)
u_bcs = FieldBoundaryConditions(top=no_slip, bottom=no_slip, north=no_slip, south=no_slip)
v_bcs = FieldBoundaryConditions(top=no_slip, bottom=no_slip)
w_bcs = FieldBoundaryConditions(north=no_slip, south=no_slip)
boundary_conditions = (u=u_bcs, v=v_bcs, w=w_bcs)

forcing = (; u=1)
closure = ScalarDiffusivity(ν=1/3300)
model = NonhydrostaticModel(; grid, closure, boundary_conditions, forcing)

simulation = Simulation(model, Δt=1e-3, stop_iteration=100)

run!(simulation)

u, v, w = model.velocities

U = Field(Average(u, dims=1))
u′ = u - U

e_op = @at (Center, Center, Center) (u′^2 + v^2 + w^2) / 2
e = Field(e_op)
compute!(e)

Ny = size(grid, 2)

using CairoMakie
heatmap(view(e, :, :, 128))


