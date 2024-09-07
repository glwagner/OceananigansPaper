using Oceananigans

grid = RectilinearGrid(topology = (Bounded, Flat, Flat),
                       size = 4,
                       x = [0, 0.1, 0.3, 0.6, 1])

c = Field{Center, Center, Center}(grid)
u = Field{Face, Center, Center}(grid)

xc = xnodes(c)
xu = xnodes(u)

scatter!(ax, xc, 0 * xc, marker=:circle, markersize=10, label="Cell centers")
scatter!(ax, xu, 0 * xu, marker=:vline, markersize=20, label="Cell interfaces")

ylims!(ax, -1, 1)
xlims!(ax, -0.1, 1.1)
hideydecorations!(ax)
hidexdecorations!(ax, ticklabels=false, label=false)
hidespines!(ax)

Legend(fig[0, 1], ax, nbanks=2, framevisible=false)

current_figure()
                       
