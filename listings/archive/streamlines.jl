using Oceananigans

function streamlines((u, v, w), positions, stop_time, sample_frequency=3)
    grid = u.grid

    x, y, z = positions
    particles = LagrangianParticles(; x, y, z)
    velocities = PrescribedVelocityFields(; u, v, w)

    model = HydrostaticFreeSurfaceModel(; grid, velocities, particles)

    umax = maximum(abs, u)
    vmax = maximum(abs, v)
    wmax = maximum(abs, w)
    Umax = max(umax, vmax, wmax)
    Δx = minimum_xspacing(grid)
    Δy = minimum_yspacing(grid)
    Δz = minimum_zspacing(grid)
    Δmin = min(Δx, Δy, Δz)

    Δt = 0.1 * Δmin / Umax
    simulation = Simulation(model; Δt, stop_time)

    # Initial position
    samples = [[(x[p], y[p], z[p]) for p = 1:length(x)]]

    function collect_sample!(sim)
        if iteration(sim) == 0
            return nothing
        end

        x = sim.model.particles.properties.x
        y = sim.model.particles.properties.y
        z = sim.model.particles.properties.z
        new_sample = [(x[p], y[p], z[p]) for p = 1:length(x)]
        push!(samples, new_sample)

        return nothing
    end

    add_callback!(simulation, collect_sample!, IterationInterval(sample_frequency))

    run!(simulation)

    # Rewrite samples[t][p] as samples[p].x
    Nt = length(samples)
    Np = length(x)

    p = 1
    xp = [samples[n][p][1] for n = 1:Nt]
    yp = [samples[n][p][2] for n = 1:Nt]
    zp = [samples[n][p][3] for n = 1:Nt]
    reshaped_samples = [(; x=xp, y=yp, z=zp)]

    for p = 2:Np
        xp = [samples[n][p][1] for n = 1:Nt]
        yp = [samples[n][p][2] for n = 1:Nt]
        zp = [samples[n][p][3] for n = 1:Nt]
        sample = (; x=xp, y=yp, z=zp)
        push!(reshaped_samples, sample)
    end

    return reshaped_samples
end


