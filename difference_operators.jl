δxᶠᶜᶜ(i, j, k, grid, c) = c[i, j, k] - c[i-1, j, k]
δxᶜᶜᶜ(i, j, k, grid, c) = c[i-1, j, k] - c[i, j, k]

ℑxᶠᶜᶜ(i, j, k, grid, c) = (c[i, j, k] + c[i-1, j, k]) / 2
ℑxᶜᶜᶜ(i, j, k, grid, u) = (u[i-1, j, k] + u[i, j, k]) / 2

δxᶠᶜᶜ(i, j, k, grid, f::Function, args...) = f(i, j, k, grid, args...) - f(i-1, j, k, grid, args...)
δxᶜᶜᶜ(i, j, k, grid, f::Function, args...) = f(i-1, j, k, grid, args...) - f(i, j, k, grid, args...)

δ²xᶜᶜᶜ(i, j, k, grid, f::Function, args...) = δxᶜᶜᶜ(i, j, k, grid, δxᶠᶜᶜ, f, args...)
