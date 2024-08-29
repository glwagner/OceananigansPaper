using Oceananigans

grid = RectilinearGrid(size = (10, 10),
                       x = (0, 1),
                       z = (0, 1),
                       topology = (Periodic, Flat, Bounded))

mask = GaussianMask{:x}(center=1, width=0.1)
sponge = Relaxation(rate=1; mask)

model = HydrostaticFreeSurfaceModel(; grid, forcing = (; u=sponge))
                                    
simulation = Simulation(model; Δt=1, stop_iteration=1)
run!(simulation)

