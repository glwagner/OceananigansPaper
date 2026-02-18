# Reviewer 5

> This paper presents an overview of the programming interface and underlying equations of the
> Oceananigans ocean-modelling framework. The overview is primarily provided through a series
> of worked examples, from "hello ocean" through to a global circulation model and even a full
> ocean-ice coupled model. For each case, the governing equations are described, their
> utilisation through Oceananigans (or ClimaOcean) is given in listing form and a sample of the
> output is visualised.

## Key Concerns

### 1. Paper aim unclear

> It is not immediately clear what the paper is aiming to communicate: whether it is primarily
> about the programming interface or about the validation of the implementation. My understanding
> is that it is about the former since there is no comparison of any of the generated results
> with those from other models.

If the paper is about the interface, then what is the purpose of the mathematical background?
If the math is sufficiently novel, its efficacy needs to be evidenced through comparison with
existing implementations. If it is not novel, why include it?

### 2. How does Oceananigans reduce complexity?

> A lot of the presented usability benefits really seem to derive from the fact that it provides
> another layer of indirection. Much of the complexity of ocean models is down to the fact that
> they necessarily have a very large number of parameters -- they are scientific instruments,
> not black boxes.

The paper would be improved by describing *how* Oceananigans reduces or helps with this
complexity. Concrete example: the forcing data lines

```julia
backend = ClimaOcean.JRA55NetCDFBackend(41)
atmosphere = ClimaOcean.JRA55_prescribed_atmosphere(arch; backend)
```

are very nice, but the paper should explain what these lines do (file download, caching, etc.)
and what a user would do if their forcing data isn't available via this mechanism.

### 3. Arithmetic and discrete calculus abstractions

> The power and expressiveness of the abstractions for arithmetic and discrete calculus could
> be very valuable in reducing developer effort while producing more reliable code.

Describe in more detail how this functionality is used.

### 4. Interface extensibility

> Ocean models are often used to do new science, so users typically want to modify them.
> Although the presented interface is nice, it would be interesting to see more description
> of what they would have to do when the interface needs to be extended.

### 5. Performance claims

> While the SYPD figures quoted for certain configurations are impressive, unqualified
> statements about Oceananigans being the "fastest ocean modelling engine" are not helpful.
> What is required is a rigorous, like-for-like comparison, ensuring the simulations are of
> comparable physical fidelity. This probably requires a wider community effort (including
> agreement on suitable benchmarks).

### 6. Tripolar grid north-fold boundary condition

> In 3.2.2, one of the examples uses a tripolar grid. How is the north-fold boundary condition
> handled? This is not covered in Silvestri (2025).
