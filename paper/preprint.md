# High-level, high-resolution ocean modeling at all scales with Oceananigans

Gregory L. Wagner<sup>1,2</sup>, Simone Silvestri<sup>1</sup>, Navid C. Constantinou<sup>3,4</sup>,
Ali Ramadhan<sup>5</sup>, Jean-Michel Campin<sup>1</sup>, Chris Hill<sup>1</sup>, Tomas Chor<sup>6</sup>,
Jago Strong-Wright<sup>7</sup>, Xin Kai Lee<sup>1</sup>, Francis Poulin<sup>8</sup>, Andre Souza<sup>1</sup>,
Keaton J. Burns<sup>1,9</sup>, Siddhartha Bishnu<sup>1,7</sup>, John Marshall<sup>1</sup>, and Raffaele Ferrari<sup>1</sup>

<sup>1</sup>Massachusetts Institute of Technology, Cambridge, MA, USA
<sup>2</sup>Aeolus Labs, San Francisco, CA, USA
<sup>3</sup>University of Melbourne, Parkville, VIC, Australia
<sup>4</sup>Australian Research Council Center of Excellence for the Weather of the 21st Century, Australia
<sup>5</sup>atdepth MRV, Cambridge, MA, USA
<sup>6</sup>University of Maryland, College Park, MD, USA
<sup>7</sup>University of Cambridge, Cambridge, United Kingdom
<sup>8</sup>University of Waterloo, Waterloo, ON, Canada
<sup>9</sup>Flatiron Institute, New York, NY, USA

## Key Points

- Oceananigans implements a programmable user interface for simulating oceanic motion at all scales.
- A high-level interface, basic numerics, and GPU-enabled high-resolution yields accessible code capable of high-fidelity simulations.
- Combining accessibility with state-of-the-art performance accelerates model development, and thus the progress of Earth system science.

## Abstract

We describe the user interface, governing equations, and numerical methods underpinning the community ocean modeling software called "Oceananigans". Oceananigans development has been lead by the Climate Modeling Alliance to build a trainable climate model with quantifiable uncertainty. Oceananigans is written in the Julia programming language, which, like similar recent efforts based on modern programming languages, distinguishes it from usual software based on Fortran. Oceananigans can efficiently simulate all scales of ocean motion, ranging from millimeter-scale turbulence in a small box to planetary-scale ocean circulation. Oceananigans design combines (i) a basic structured finite volume algorithm (ii) optimized for high-resolution simulations on GPUs which is (iii) exposed behind a high-level, programmable user interface. This design negotiates a dual mandate for highest-possible performance (to support state-of-the-art applications) and enhanced accessibility (to facilitate adoption and development). The dual mandate aims ultimately to accelerate the progress of Earth system science. Achieving this aim, however, requires a substantial and sustained increase in the collective effort of Oceananigans development.

## Plain Language Summary

This paper describes Oceananigans, a tool for simulating ocean currents and fluid motion written in the Julia programming language. Using a relatively new programming language separates Oceananigans from usual software written in Fortran and places it alongside a handful of recent efforts to modernize Earth system modeling software. Oceananigans can simulate ocean motion on a wide range of scales from millimeters to planetary-scale. Oceananigans is also the fastest ocean physics simulator to date, because it was written from scratch for graphics processing units (GPUs). We argue that the Oceananigans modeling strategy, which combines basic numerics on GPUs with a powerful user interface, can accelerate the pace of model development and therefore progress in Earth system science.

## 1 Introduction

Only numerical models can synthesize the vast accumulation of ocean and Earth system knowledge. As a result, the capabilities of modeling software rate-limit progress in Earth system science. Since the first general circulation models ran on primitive computers (Phillips, 1956; Bryan, 1969), advances in hardware, numerical methods, and the approximate parameterization of otherwise unresolved processes have improved the fidelity of ocean simulations (Griffies et al., 2015). But the gap between potential and practice in ocean modeling is widening because advances in software and hardware technology have outpaced model development. Today, most ocean modeling software (i) does not run on the fastest computers (which are based on GPUs), (ii) relies on outdated user interfaces that are less efficient than modern approaches and unfamiliar to a new generation of programmers, and (iii) can tackle only a limited subset of ocean modeling problems.

This paper describes new ocean modeling software written in the Julia programming language (Bezanson et al., 2017) called Oceananigans. Oceananigans is being developed by the Climate Modeling Alliance and external collaborators as part of a larger effort to develop a climate model automatically-calibrated to observations and high resolution simulations, and with quantified uncertainty (Schneider et al., 2017). This goal first and foremost drives Oceananigans design: fast enough for ensemble-based calibration (Silvestri et al., 2025), furnished with trainable parameterizations (Wagner, Hillier, et al., 2025), and flexible enough to support a hierarchical modeling approach (Held, 2005) that requires simulations ranging from small-scale days-long large eddy simulations to global climate simulations over decades.

Oceananigans development also strives to close the gap between potential and practice by prototyping a framework that can accelerate the process of model development. To this end, we describe a numerical approach to ocean modeling that: (i) can model all scales of oceanic motion at high resolution, (ii) enables straightforward algorithmic and performance optimization tailored to modern hardware, and (iii) facilitates a user interface design that takes maximum advantage of the use of a modern programming languages. We hope that faster development of models like Oceananigans will ultimately, through a longer process of collective effort, accelerate progress in ocean and climate science. Yet our progress is incomplete, and more model development, further improvements to our modeling framework, and continued growth of the Oceananigans community will be required to realize this objective.

### 1.1 Modeling motions from millimeters to millennia

The evolution of ocean circulation over millennia is controlled by turbulent mixing with scales that range down to millimeters. Two distinct systems have evolved to model this huge range of oceanic motion: "GCMs" (general circulation models) for hydrostatic regional-to-global scale simulations, and simpler software for nonhydrostatic large eddy simulations (LESs) with meter-scale resolution that are high-fidelity but limited in duration and extent. Compared to LESs, GCMs usually invoke more elaborate numerical methods and parameterizations to cope with the global ocean's complex geometry and the more significant impacts of unresolved subgrid processes.

Oceananigans began as software for LES (Ramadhan et al., 2020), by refining an approach for hybrid hydrostatic/nonhydrostatic dynamical cores pioneered by MITgcm (Marshall et al., 1997) for GPUs. Our nonhydrostatic LES algorithm was then adapted and optimized for a hydrostatic GCM (Silvestri et al., 2025), yielding significant acceleration over existing CPU-based codes (Silvestri et al., 2023). At the same time, we developed LES-inspired, minimally dissipative numerical methods for turbulence-resolving simulations on finite volume C-grids (Silvestri et al., 2024) with up to 2x higher effective resolution than existing advection schemes, and which automatically adapt local flow conditions and grid spacing.

The result is an efficient modeling system suited to a brute force, resolution-focused approach to accuracy for all scales of oceanic motion. For example, Silvestri et al. (2023) and Silvestri et al. (2025) report near-global, ocean-only simulations achieving 10 SYPD on 64 A100s at 8 km resolution and 1 SYPD on 512 A100s at 2 km resolution (using a now-outdated version of Oceananigans). Section 3.2.3 shows results from a global coupled ocean and sea ice simulation with 16 km resolution which, with the current version of Oceananigans, achieves 1.5 simulated years per day on a single Nvidia H100 GPU (representing efficiency improvements to Oceananigans as well as the superior performance of H100s versus A100s). In comparison, the 25-km-resolution OM4 configuration reported in Adcroft et al. (2019) would achieve the same throughput on 22 nodes of their Cray XC40 machine as our 16-km-resolution Oceananigans simulation on a single H100. If we estimate that a single XC40 node is equivalent to 1--4 H100s, these numbers imply that Oceananigans provides an efficiency gain of 20--80x over the state of the art. Silvestri et al. (2025) argue that Oceananigans efficiency enables routine climate simulations that resolve oceanic mesoscale turbulence rather than parameterizing it, eliminating a major bias affecting the accuracy of climate projections.

Such an "LES the ocean" strategy has the advantage of simplicity over strategies that require sophisticated models for explicit dissipation, generalized vertical coordinates (Shchepetkin & McWilliams, 2005; Leclair & Madec, 2011; Petersen et al., 2015), Lagrangian vertical advection (Halliwell, 2004; Griffies et al., 2020), or unstructured horizontal grids (Ringler et al., 2013; Danilov et al., 2017; Korn et al., 2022). Moreover, high resolution yields a plethora of additional improvements in simulation fidelity (Chassignet & Xu, 2017; Kiss et al., 2020; Chassignet & Xu, 2021). A major open question regarding the choice of vertical coordinate, however, is whether higher resolution will also reduce the spurious numerical mixing that pollutes the fidelity of lower-resolution simulations with simple vertical coordinates (e.g. Griffies et al., 2000).

### 1.2 How modern software can accelerate progress in ocean and climate science

In addition to the goal of providing a performant nonhydrostatic-to-hydrostatic ocean modeling system, Oceananigans also attempts to accelerate model development by using a modern programming language. In traditional paradigms based on Fortran, for example, new parameterizations and numerical methods are typically prototyped and validated in a "productivity" language such as MATLAB or Python, prior to implementation in a production context. With the Julia programming language, on the other hand, Oceananigans provides a framework wherein new parameterizations (Wagner, Hillier, et al., 2025), new numerical methods (Silvestri et al., 2024), and new algorithms for performance optimization (Silvestri et al., 2025) may be both prototyped and refined for production without re-implementation. We note that this kind of framework may also be implemented with a Python-based domain-specific languages, such as JAX (Bradbury et al., 2018). Modern programming languages increase the productivity of users in addition to model developers.

The need for new frameworks to accelerate progress in Earth system modeling is evidenced by the recent proliferation of similar efforts: for example, regional-mom6 (Barnes et al., 2024) and CROCO-tools (Jullien et al., 2025) provide Python software that automates the configuration of existing models (Adcroft et al., 2019; Auclair et al., 2025) by user-scientists. Such tools provide advantages over new modeling systems by leveraging mature modeling systems with massive user communities and decades of development history.

New modeling systems are also being developed in modern programming languages. ClimaAtmos (Yatunin et al., 2025) and ClimaLand, for example, provide software for nonhydrostatic atmosphere simulations and land surface modeling implemented in pure Julia. Veros (hydrostatic ocean simulations, Hafner et al., 2021; Mrozowska et al., 2025) and NeuralGCM (hydrostatic atmosphere simulations, Kochkov et al., 2024; Yuval et al., 2024) are written in Python using JAX (Bradbury et al., 2018). JAX-based software benefits from the ubiquity and depth of the Python ecosystem and offers portable performance between CPUs, GPUs, and TPUs. Most notably, JAX supports automatic differentiation, enabling gradient-based parameter estimation or training of machine learning (ML) components (Kochkov et al., 2021, 2024) and adjoint-based data assimilation (Solvik et al., 2025). While Oceananigans achieves performance portability via KernelAbstractions (Churavy, 2024), and differentiability is being developed via Enzyme and Reactant (Moses et al., 2021), neither tool is as mature as JAX.

In outward appearances -- its user interface and primary novelty -- Oceananigans implements a programmable, script-first user interface for configuring, running, monitoring, and analyzing simulations. Oceananigans' user interface design was inspired primarily by Dedalus (Burns et al., 2020), Python software that solves PDEs parsed from strings using global spectral methods. Programmable interfaces are also implemented by SpeedyWeather (Klower et al., 2024), Julia-based software for hydrostatic atmosphere simulations using global spectral methods, and Thetis (Karna et al., 2018), which leverages Firedrake (Rathgeber et al., 2016) to implement the hydrostatic primitive equations for ocean modeling using discontinuous Galerkin numerical methods.

### 1.3 Why programmable interfaces matter

In 1984, Cox published the first description of generalizable ocean modeling software (Cox, 1984; Griffies et al., 2015). The "Cox model" is written in FORTRAN 77 and features a multi-step user interface for building new models: first, source code modifications are written to determine, for example, domain geometry and boundary conditions, emplaced into the "base code", and compiled. Next, in a second step, a text-based namelist file is used to determine parameters like the stop iteration, mixing coefficients, and solver convergence criteria. Cox (1984) provided three examples to illustrate the user interface, providing both source code and namelists for each example.

With more than forty years of progress in software engineering, numerical methods, and parameterization of unresolved processes, and more than a billion times more computational power, most of today's ocean models bear little resemblance to the Cox model -- except, perhaps, for their user interfaces. Most ocean models still invoke the non-programmable, namelist-based paradigm described by the Cox model documentation -- even modeling systems written in modern programming languages. Model workflows often involve multiple steps to generate input data, configure a set of namelist files, modify source code to change the model equations in ways not accessible through a change of parameters, and finally to compile and run and model software.

A central thesis of this paper is that improvements to user interfaces to Earth system modeling software are essential for accelerating progress in Earth system science. In particular, a programmable user interface can provide a seamless one-step workflow for numerical experiments including setup, execution, analysis, and visualization with a single script. Programmable interfaces written in scripting languages like Python and Julia are the engine of progress (Perez & Granger, 2007) in fields ranging from visualization (e.g. matplotlib, Hunter, 2007), to machine learning (e.g. pytorch, Paszke et al., 2019), to physics (e.g. dedalus, fenics, or firedrake, Alnaes et al., 2015; Rathgeber et al., 2016; Burns et al., 2020). Programmable interfaces facilitate fast prototyping, collaboration through code sharing, and reproducible simulations with a small number of files.

Oceananigans implements a programmable, library-style interface to Earth system modeling software written in the Julia programming language. We emphasize that programmable interfaces can be implemented in any scripting language (all of the cited examples above are based on Python), and we do not argue that Julia is strictly superior to other approaches. That said, it bears mentioning some of the advantages and disadvantages of using Julia compared to approaches based on Python/JAX. Python's main advantage is its status as the lingua franca of programming languages: Python has a much larger open source community, and greater resources devoted to its support and development. Some of these advantages are mitigated by the seamless interoperability between Julia and Python, the advent of powerful AI-based coding tools, and the similarity between Julia and Python syntax. Both Julia and Python enable interactivity, extensibility, automatic installation on any system, and portability to laptops and GPUs (Besard et al., 2018; Bradbury et al., 2018; Churavy, 2024). The main advantage of Julia is that high-performance code may be developed without relying on a domain-specific language such as JAX. This yields several unique capabilities: first, a wider range of strategies can be easily deployed for performance optimization, or for high performance implementations of unusual ocean-modeling-specific algorithms (Besard et al., 2018). Second, scripted user code that implements custom forcing and boundary conditions is easily embedded within a simulation that runs on GPUs (Besard et al., 2018).

### 1.4 Outline of this paper

This paper proceeds in section 2 by illustrating the basic form of Oceananigans' programmable interface with two classroom examples: two-dimensional turbulence, and a passive tracer advected by two-dimensional turbulence and forced by user-specified forcing. Our goal is to demonstrate how Oceananigans' user interface both simplifies basic simulations while enabling complex, creative science. We do not attempt to document the specifics of the user interface in detail or to provide a comprehensive description of all features, however: for that we refer the reader to Oceananigans documentation.

Section 3 continues by writing down the Boussinesq governing equations that underpin Oceananigans' nonhydrostatic and hydrostatic models. We then sketch out Oceananigans capabilities with a series of examples that progress from basic direct numerical simulations of cabbeling, to realistic tidally-forced nonhydrostatic large eddy simulations over a headland, to a 1/6th degree hydrostatic, eddying, global ocean simulation. To summarize, Oceananigans supports both nonhydrostatic and hydrostatic simulations on rectilinear and curvilinear grids with or without bathymetry. A suite of pressure solvers, WENO-based advection schemes (C. Shu, 1997), Laplacian diffusivity, and subgrid closures including constant and dynamic Smagorinsky (Smagorinsky, 1963; Lilly, 1983; Bou-Zeid et al., 2005) and Anisotropic Minimum Dissipation (Rozema et al., 2015; Vreugdenhil & Taylor, 2018) support nonhydrostatic direct and large eddy simulations of meter-scale phenomena in periodic, closed, or open domains. Hydrostatic regional to global ocean simulations with an implicit or split-explicit free surface formulation and z or z-star vertical coordinate are supported on rectilinear, latitude-longitude, tripolar (Murray, 1996), and cubed sphere (Adcroft et al., 2004) grids with second-order and WENO-based vector invariant schemes (Silvestri et al., 2024), the Gent--McWilliams parameterization (Gent & Mcwilliams, 1990), horizontal biharmonic diffusivity, and Ri-based, one-equation (Wagner, Hillier, et al., 2025), and two-equation vertical mixing schemes (Umlauf & Burchard, 2005). Linear, quadratic (Roquet, Madec, Brodeau, & Nycander, 2015), and polynomial (Roquet, Madec, McDougall, & Barker, 2015) equations of state are supported (with arbitrary gravitational direction for nonhydrostatic simulations) as well as traditional, non-traditional, spherical, and beta-plane Coriolis forces. Lagrangian particles, Stokes drift, biogeochemistry (Strong-Wright et al., 2023) are supported. A second package called "ClimaOcean" provides a coupled modeling framework for integrating Oceananigans ocean models together with sea ice simulations and either prescribed or prognostic atmospheres, using fluxes computed according to Monin-Obukhov similarity theory (for example, Edson et al., 2014). ClimaOcean also implements utilities for automatic, script-based downloading of bathymetric data, ocean and sea ice physical and biogeochemical reanalysis products, and atmospheric reanalysis. (A name change to ClimaOcean is underway, as the package is evolving into a comprehensive tool for coupled Earth system modeling and is no longer restricted to ocean modeling.)

We conclude in section 4 by describing outstanding problems and questions, outlining future development work, and anticipating the next major innovations in ocean modeling that will someday render the present work obsolete. Ongoing work includes coupling with atmosphere models (Yatunin et al., 2025; Klower et al., 2024), and the development of differentiable workflows with Enzyme and Reactant (Moses et al., 2021) that invoke the integration of Oceananigans models. Important future work includes the continued implementation, development, and calibration of more theory-based and ML-based parameterizations for vertical mixing (Harcourt, 2015; Reichl & Li, 2019; Legay et al., 2025; Wagner & Constantinou, 2025), mesoscale turbulence (Mak et al., 2018; Jansen et al., 2019), bathymetry-mesoscale interaction, ice-ocean coupling, and air-sea coupling (Pelletier et al., 2021).

## 2 Oceananigans, the library

Oceananigans is fundamentally a library of tools for building models by writing programs called "scripts". This design departs from typical monolithic interfaces that ingest lists of flags and parameters from non-executable text files. By blending mathematical symbols with verbose natural language names, Oceananigans syntax tries to enable evocative scripting that approaches the effectiveness of writing for communicating computational science.

### 2.1 Hello, ocean

Learning Oceananigans starts with running simple simulations. Our first example in listing 1 sets up, runs, and visualizes a simulation of two-dimensional turbulence. The 21 lines of listing 1 illustrate one of Oceananigans' main achievements: a numerical experiment may be completely described by a single script. To execute the code in listing 1, we need to copy into a file (call this, for example, hello_ocean.jl) and executed by typing `julia hello_ocean.jl` at a terminal.

**Listing 1:** A Julia script that uses Oceananigans and the Julia plotting library CairoMakie to set up, run, and visualize a simulation of two-dimensional turbulence on a Graphics Processing Unit (GPU). The initial velocity field, defined on lines 11-12, consists of random numbers uniformly-distributed between -1 and 1. The vorticity $\zeta = \partial_x v - \partial_y u$ is defined on line 18. The solution is visualized in figure 1.

```julia
using Oceananigans, CUDA # using CUDA allows us to use an Nvidia GPU

# The third dimension is "flattened" to reduce the domain from three to two dimensions.
topology = (Periodic, Periodic, Flat)
architecture = GPU() # CPU() works just fine too for this small example.
x = y = (0, 2pi)
grid = RectilinearGrid(architecture; size=(256, 256), x, y, topology)

model = NonhydrostaticModel(; grid, advection=WENO(order=9))

epsilon(x, y) = 2rand() - 1 # Uniformly-distributed random numbers between [-1, 1).
set!(model, u=epsilon, v=epsilon)

simulation = Simulation(model; dt=0.01, stop_time=10)
run!(simulation)

u, v, w = model.velocities
zeta = dx(v) - dy(u)

using CairoMakie
heatmap(zeta, colormap=:balance, axis=(; aspect=1))
```

Oceananigans scripts organize into four sections. The first three define the grid, model, and simulation, and conclude with execution of the simulation. The fourth section, often implemented separately for complex or expensive simulations, performs post-processing and analysis. In listing 1, the grid defined on lines 4--7 determines the problem geometry, spatial resolution, and machine architecture. To use a CPU instead of a GPU, one writes `CPU()` in place of `GPU()` on line 5: no other changes to the script are required.

Lines 9--12 define the model, which solves the Navier--Stokes equations in two dimensions with a 9th-order Weighted, Essentially Non-Oscillatory (WENO) advection scheme (see section Appendix B and Silvestri et al. (2024) for more information about WENO). The velocity components u, v are initialized with uniformly distributed random numbers within [-1, 1). The model definition can also encompass forcing, boundary conditions, and the specification of additional terms in the momentum and tracer equations such as Coriolis forces or turbulence closures.

Line 14 builds a simulation with a time-step $\Delta t = 0.01$ which will run until $t = 10$ (time is non-dimensional via user input in this case). Lines 17-18 analyze the final state of the simulation by computing vorticity, illustrating Oceananigans' toolbox for building expression trees of discrete calculus and arithmetic operations. The same tools may be used to define online diagnostics to be periodically computed and saved to disk while the simulation runs. Line 21 concludes the numerical experiment with a visualization, which is shown in figure 1.

### 2.2 Incorporating user code

With a programmable interface and aided by Julia's just-in-time compilation, user functions specifying domain geometry, forcing, boundary conditions, and initial conditions can be incorporated directly into models without a separate programming environment. To illustrate function-based forcing, we modify listing 1 with code that adds a passive tracer which is forced by a moving source that depends on x, y, t. A visualization of the vorticity and tracer field generated by listings 1 and 2 are shown in figure 1.

**Listing 2:** Implementation of a moving source of passive tracer with a function in a two-dimensional turbulence simulation. These lines of code replace the model definition on line 9 in listing 1.

```julia
function circling_source(x, y, t)
    delta, omega, r = 0.1, 2pi/3, 2
    dx = x + r * cos(omega * t)
    dy = y + r * sin(omega * t)
    return exp(-(dx^2 + dy^2) / 2delta^2)
end

forcing = (; c = circling_source)
model = NonhydrostaticModel(; grid, advection=WENO(order=9), tracers=:c, forcing)
```

[Figure 1: Vorticity after t = 10 (left) and a passive tracer injected by a moving source at t = 2.5 (right) in a simulation of two-dimensional turbulence using an implicitly-dissipative advection scheme.]

General tasks may be inserted into the time stepping loop by modifying the simulation. This supports things as mundane as printing a summary of the current model status or writing output, to more exotic tasks like nudging state variables or updating a diffusion coefficient based on an externally-implemented model.

### 2.3 Abstractions for arithmetic and discrete calculus

Abstractions representing unary, binary, and calculus operators produce a system for building "lazy" expression trees that are evaluated only upon request (for example, if their output should be periodically saved to disk during a simulation). Example calculations representing vorticity, $\zeta = \partial_x v - \partial_y u$, speed $s = \sqrt{u^2 + v^2}$, and the x-integral of enstrophy $Z = \int_0^{2\pi} \zeta^2 \, dx$ are shown in listing 3.

**Listing 3:** Illustration of abstractions for expression trees and reductions called "AbstractOperations" in Oceananigans. AbstractOperations are "lazy", in the sense that they represent potential computations but do not instantiate the computation directly. One of the main use cases for AbstractOperations is to support diagnostic computations performed repeatedly throughout the course of a simulation.

```julia
u, v, w = model.velocities

# "Lazy" expression trees and reductions representing computations:
zeta = dx(v) - dy(u)
s = sqrt(u^2 + v^2)
Z = Integral(zeta^2, dims=1)

# Building and computing a Field that instantiates an AbstractOperation:
zeta_field = Field(dx(v) - dy(u))
compute!(zeta_field)
```

## 3 Governing equations, parameterizations, and illustrative examples

Oceananigans implements two "models" for ocean-flavored fluid dynamics: the HydrostaticFreeSurfaceModel, and the NonhydrostaticModel. Each represents a template for equations that govern the evolution of momentum and tracers. Both models are incompressible and make the Boussinesq approximation, which means that the density of the modeled fluid is decomposed into a constant reference $\rho_0$ and a small dynamic perturbation $\rho'$,

$$\rho(\mathbf{x}, t) = \rho_0 + \rho'(\mathbf{x}, t) \quad \text{where} \quad \rho' \ll \rho_0 \,, \tag{1}$$

and $\mathbf{x} = (x, y, z)$ is position and $t$ is time.

The relative smallness of $\rho'$ reduces conservation of mass to a statement of incompressibility called the continuity equation,

$$\nabla \cdot \mathbf{u} = 0 \,, \tag{2}$$

where

$$\mathbf{u} \stackrel{\text{def}}{=} u \, \hat{x} + v \, \hat{y} + w \, \hat{z} \,, \tag{3}$$

is the three-dimensional velocity field. Within the Boussinesq approximation, the momentum $\rho_0 \mathbf{u}$ varies only with the velocity $\mathbf{u}$. The effect of density variations is encapsulated by a buoyant acceleration,

$$b \stackrel{\text{def}}{=} -\frac{g \rho'}{\rho_0} \,, \tag{4}$$

where $g$ is gravitational acceleration. The "buoyancy" $b$ acts in the direction of gravity.

The total dynamic pressure $P$ is decomposed into

$$P = \rho_0 g z + \rho_0 p(\mathbf{x}, t) \,, \tag{5}$$

where $z$ is height, $\rho_0 g z$ is the static contribution to pressure that opposes the gravitational force associated with the reference density $\rho_0$, and $\rho_0 p$ represents the dynamic anomaly. $p$ is called the kinematic pressure.

### 3.1 The NonhydrostaticModel

The NonhydrostaticModel represents the Boussinesq equations formulated without making the hydrostatic approximation typical to general circulation models. The NonhydrostaticModel has a three-dimensional prognostic velocity field.

#### 3.1.1 The NonhydrostaticModel momentum equation

The NonhydrostaticModel's momentum equation incorporates advection by a background velocity field, Coriolis forces, surface wave effects via the Craik-Leibovich asymptotic model (Craik & Leibovich, 1976; Huang, 1979), a buoyancy term allowed to be a nonlinear function of tracers and depth, a stress divergence derived from molecular friction or a turbulence closure, and a user-defined forcing term. Using the Boussinesq approximation in (1) and the pressure decomposition in (5), the generic form of NonhydrostaticModel's momentum equation is

<!-- [Equation 6: NonhydrostaticModel momentum equation] -->
$$\partial_t \mathbf{u} = -\nabla p \underbrace{- (\mathbf{u} \cdot \nabla)\mathbf{u} - (\mathbf{u}_b \cdot \nabla)\mathbf{u} - (\mathbf{u} \cdot \nabla)\mathbf{u}_b}_{\text{advection}} \underbrace{- \mathbf{f} \times \mathbf{u}}_{\text{Coriolis}} + \underbrace{(\nabla \times \mathbf{u}_s) \times \mathbf{u} + \partial_t \mathbf{u}_s}_{\text{Stokes forcing}} \underbrace{- b\,\hat{g}}_{\text{buoyancy}} \underbrace{- \nabla \cdot \boldsymbol{\tau}}_{\text{closure}} + \underbrace{\mathbf{F}_u}_{\text{forcing}} \,, \tag{6}$$

where $\mathbf{u}_b$ is a prescribed and divergence-free "background" velocity field, $p$ is the kinematic pressure, $\mathbf{f}$ is the background vorticity associated with a rotating frame of reference, $\mathbf{u}_s$ is the Stokes drift profile associated with a prescribed surface wave field, $b$ is buoyancy, $\hat{g}$ is the gravitational unit vector (usually pointing downwards, that is, $\hat{g} = -\hat{z}$), $\boldsymbol{\tau}$ is the stress tensor associated with molecular viscous or subgrid turbulent momentum transport, and $\mathbf{F}_u$ is a body force.

To integrate equation (6) while enforcing (2), we use a pressure correction method that requires solving a three-dimensional Poisson equation to find $p$, which can be derived from $\nabla \cdot$ (6). This Poisson equation is often a computational bottleneck in curvilinear or irregular domains, and its elimination is the main motivation for making the hydrostatic approximation when formulating the HydrostaticFreeSurfaceModel, as described in section 3.2. For rectilinear grids, we solve the Poisson equation using a direct FFT-based or mixed FFT-tridiagonal solver (Schumann & Sweet, 1988), providing substantial acceleration over MITgcm's conjugate gradient pressure solver (Marshall et al., 1997). In irregular domains, we either use a masking method that permits an approximate solution of the pressure Poisson equation with the FFT-based method, or an iterative conjugate gradient solver that leverages the FFT-based solver as a preconditioner. The pressure correction scheme is described further in appendix A2.

Using (2), advection in the NonhydrostaticModel may formulated in "flux form",

$$\text{advection} = u_j \partial_j u_i + u_{bj} \partial_j u_i + u_j \partial_j u_{bi} = \partial_j \left[ (u_j + u_{bj}) u_i + u_j u_{bi} \right] \,, \tag{7}$$

where, we have used indicial notation and for example, the i-th component of the advection term is $[(\mathbf{u} \cdot \nabla)\mathbf{u}]_i = u_j \partial_j u_i$. (See the text surrounding equations (18)--(20) for a discussion of advection term formulation in the HydrostaticFreeSurface model.)

The formulation of the Stokes drift terms means that $\mathbf{u}$ is the Lagrangian-mean velocity when Stokes drift effects are included (see, for example, Wagner et al., 2021). With a Lagrangian-mean formulation, equations (2) and (6) are consistent only when $\mathbf{u}_s$ is non-divergent -- or equivalently, when $\mathbf{u}_s$ is obtained by projecting the divergence out of the usual Stokes drift (Vanneste & Young, 2022). As discussed by Wagner et al. (2021), the Lagrangian-mean formulation of (6) means that closures for LES strictly destroy kinetic energy, avoiding the inconsistency between resolved and subgrid fluxes affecting typical LES formulated in terms of the Eulerian-mean velocity (see also Pearson, 2018; Wagner & Constantinou, 2025).

The labeled terms in (6) are controlled by arguments to NonhydrostaticModel invoked in both of listings 1 and 2. For example, "advection" chooses a numerical scheme to approximate the advection term in (6) and (7). As another example, we consider configuring the closure term in (6) to represent (i) molecular diffusion by a constant-coefficient Laplacian ScalarDiffusivity, (ii) turbulent stresses approximated by the SmagorinskyLilly eddy viscosity model (Smagorinsky, 1963; Lilly, 1983) for large eddy simulation, or (iii) omitting it entirely and relying on WENO advection schemes for dissipation and stability (Pressel et al., 2017; Wagner, Hillier, et al., 2025). In these three cases, the closure flux divergence $\nabla \cdot \boldsymbol{\tau} = \partial_m \tau_{nm}$ in indicial notation becomes

<!-- [Equation 8: Closure options] -->
$$-\partial_m \tau_{nm} = \begin{cases} \partial_m (\nu \, \partial_m u_n) & \text{(ScalarDiffusivity)} \\ 0 & \text{(nothing)} \\ \partial_m \left( \underbrace{2\,C_s \Delta^2 |\Sigma|}_{\nu_e} \Sigma_{nm} \right) & \text{(SmagorinskyLilly)} \end{cases} \tag{8}$$

where $\nu$ is the Laplacian diffusion coefficient, $\Sigma_{nm} = \partial_m u_n + \partial_n u_m$ is the strain rate tensor, $|\Sigma|$ is the magnitude of the strain rate tensor, $C_s$ is the SmagorinskyLilly model constant, $\Delta$ scales with the local grid spacing, and $\nu_e$ is the eddy viscosity.

ScalarDiffusivity diffusion coefficients may also vary in time- and space or depend on model fields. Vertically-implicit time-discretization is supported with ScalarDiffusivity. Other closure options for NonhydrostaticModel include fourth-order ScalarBiharmonicDiffusivity, various Lagrangian-averaged or directionally-averaged flavors of DynamicSmagorinsky (Bou-Zeid et al., 2005), and the AnisotropicMinimumDissipation turbulence closure (Rozema et al., 2015; Abkar et al., 2016; Vreugdenhil & Taylor, 2018; Wagner et al., 2021) for large eddy simulations.

We note that large eddy simulations may be conducted solely with WENO advection schemes that dissipate grid-scale kinetic energy undergoing a forward cascade from large to small scales. However, no reliable method has yet been developed to diagnose the dissipation of kinetic energy in such simulations. This means that explicit closures must be included to diagnose kinetic energy dissipation. Anecdotally, simulations with explicit closures tend to dissipate more kinetic energy than simulations that rely purely on implicit dissipation via WENO advection (for example, Pressel et al., 2017). But more work is needed to investigate the fidelity of explicit versus implicit dissipation, and moreover to develop methods for diagnosing implicit kinetic energy and tracer variance dissipation by WENO advection.

**Listing 4:** Direct numerical simulation of flow past a cylinder at various Reynolds numbers Re. The domain is periodic in x and a sponge layer on the right side relaxes the solution to $\mathbf{u} = u_\infty \hat{x}$ with $u_\infty = 1$. The experiment can be converted to a large eddy simulation (thereby sending $Re \to \infty$) by replacing the no-slip boundary conditions with an appropriate drag model and either (i) using an appropriate turbulence closure or (ii) using the WENO(order=9) advection scheme with no turbulence closure. Visualizations of the DNS and LES cases are shown in figure 2.

```julia
r, U, Re, Ny = 1/2, 1, 1000, 2048

grid = RectilinearGrid(GPU(), size=(2Ny, Ny), x=(-3, 21), y=(-6, 6),
    topology=(Periodic, Bounded, Flat))

cylinder(x, y) = (x^2 + y^2) <= r^2
grid = ImmersedBoundaryGrid(grid, GridFittedBoundary(cylinder))

closure = ScalarDiffusivity(nu=1/Re)

no_slip = FieldBoundaryConditions(immersed=ValueBoundaryCondition(0))
boundary_conditions = (u=no_slip, v=no_slip)

# Implement a sponge layer on the right side of the domain that
# relaxes v -> 0 and u -> U over a region of thickness delta
@inline mask(x, y, delta=3, x0=21) = max(zero(x), (x - x0 + delta) / delta)
Fu = Relaxation(target=U; mask, rate=1)
Fv = Relaxation(target=0; mask, rate=1)

model = NonhydrostaticModel(; grid, closure, boundary_conditions, forcing=(u=Fu, v=Fv))
```

Listing 4 implements a direct numerical simulation of uniform flow past a cylinder with no-slip boundary conditions, a molecular ScalarDiffusivity, and a centered second-order advection scheme. Lines 6--7 embed a cylindrical mask in a RectilinearGrid using a GridFittedBoundary, which generalizes to arbitrary three-dimensional shapes. The no-slip condition is specified on lines 11--12 by invoking ValueBoundaryCondition (a synonym for "Dirichlet" boundary condition): ValueBoundaryCondition enforces a near-wall tangential velocity gradient, which in turn produces a viscous flux associated with the ScalarDiffusivity viscosity (resolving to the same formulation proposed in a seminal paper by Harlow et al., 1965). Other choices include GradientBoundaryCondition (Neumann), FluxBoundaryCondition (direct imposition of fluxes), and OpenBoundaryCondition (for non-trivial boundary-normal velocity fields).

Results obtained with listing 4 for Re = 100, Re = 1000, and a modified version of listing 4 for large eddy simulation (Re -> infinity) are visualized in figure 2. To adapt listing 4 for LES, the closure is eliminated in favor of a 9th-order WENO advection scheme, and the no-slip boundary condition is replaced with a quadratic drag boundary condition with a drag coefficient estimated from similarity theory using a constant estimated roughness length.

[Figure 2: Vorticity snapshots in simulations of flow around a cylinder. The top two panels show vorticity in direct numerical simulations (DNS) that use a molecular ScalarDiffusivity closure and Centered(order=2) advection. The bottom panel shows a large eddy simulation (LES) with no closure and a WENO(order=9) advection scheme.]

#### 3.1.2 The NonhydrostaticModel tracer conservation equation

The buoyancy term in (6) requires tracers, and can be formulated to use buoyancy itself as a tracer, or to depend on temperature $T$ and salinity $S$. For seawater, a 54-term polynomial approximation TEOS10EquationOfState (McDougall & Barker, 2011; Roquet, Madec, McDougall, & Barker, 2015) is implemented in the auxiliary package SeawaterPolynomials, along with quadratic approximations to TEOS-10 (Roquet, Madec, Brodeau, & Nycander, 2015) and a LinearEquationOfState. All tracers -- either "active" tracers required to compute the buoyancy term, as well as additional user-defined passive tracers -- obey the tracer conservation equation

<!-- [Equation 9: Tracer conservation equation] -->
$$\partial_t c = \underbrace{-\nabla \cdot [(\mathbf{u} + \mathbf{u}_b + \mathbf{u}_p)\,c + \mathbf{u}\,c_b] + c\,\nabla \cdot \mathbf{u}_p}_{\text{advection}} \underbrace{-\nabla \cdot \mathbf{J}_c}_{\text{closure}} + \underbrace{S_c}_{\text{biogeochemistry}} + \underbrace{F_c}_{\text{forcing}} \,, \tag{9}$$

where $c$ represents any tracer, $c_b$ represents a prescribed background tracer concentration for $c$, $\mathbf{u}_p$ is the velocity associated with biogeochemical transformations (and which, in general, can have non-zero divergence), $\mathbf{J}_c$ is a tracer flux associated with molecular diffusion or subgrid turbulence, $S_c$ is a source or sink term associated with biogeochemical transformations (provided, for example, by external packages like OceanBioME; Strong-Wright et al., 2023), and $F_c$ is a user-defined source or sink. The formulation of tracer advection in (9) leverages $\nabla \cdot \mathbf{u} = 0$ and $\nabla \cdot \mathbf{u}_b = 0$.

A simulation with a passive tracer having a user-defined source term is illustrated by listing 2 and figure 1. For a second example, we consider freshwater cabbeling (see, for example, Bisits et al., 2025). Cabbeling occurs when two water masses of similar density mix to form a new water mass which, due to the nonlinearity of the equation of state, is denser than either of its constituents. Freshwater, for example, is densest at 4 degrees Celsius, while 1- and 7.55-degree water are lighter with roughly the same density. We implement a direct numerical simulation (DNS) in which 7.55-degree water overlies 1-degree water, using the TEOS10EquationOfState provided by the auxiliary package SeawaterPolynomials. The script is shown in listing 5. The resulting density and temperature fields after 1 minute of simulation are shown in figure 3. Note that the TEOS10EquationOfState typically depends on both temperature and salinity tracers, but listing 5 specifies a constant salinity S = 0 and thus avoids allocating memory for or simulating salinity directly. Also note, DNS is not Oceananigans' strength due to its second-order finite volume formulation, compared to pseudospectral formulations (Lecoanet et al., 2016). (This contrasts with Oceananigans' capabilities for large eddy simulation, where simulation error is dominated by the representation of grid-scale dissipation rather than the formal accuracy of the discretization.) Future work to improve DNS with Oceananigans could investigate higher-order discretizations of viscous diffusion and molecular tracer diffusion.

**Listing 5:** Direct numerical simulation of convective turbulence driven by cabbeling between 1- and 7.55-degree freshwater. $\nu$ denotes viscosity and $\kappa$ denotes the tracer diffusivity. The diffusivity may also be set independently for each tracer.

```julia
grid = RectilinearGrid(GPU(), topology = (Bounded, Flat, Bounded),
    size = (4096, 1024), x = (0, 2), z = (-0.5, 0))

closure = ScalarDiffusivity(nu=1.15e-6, kappa=1e-7)

using SeawaterPolynomials: TEOS10EquationOfState
equation_of_state = TEOS10EquationOfState(reference_density=1000)

buoyancy = SeawaterBuoyancy(gravitational_acceleration = 9.81;
    constant_salinity = 0, # set S=0 and simulate T only
    equation_of_state)

model = NonhydrostaticModel(; grid, buoyancy, closure, tracers=:T)

Ti(x, z) = z > -0.25 ? 7.55 : 1
Xi(x, z) = 1e-2 * randn()
set!(model, T=Ti, u=Xi, v=Xi, w=Xi)
```

[Figure 3: Density and temperature at t = 1 minute in a direct numerical simulation of cabbeling in freshwater. Note that both fields span from x = 0 to x = 2 meters; only the left half of the density field and the right half of the temperature field are shown.]

We next consider a large eddy simulation of the Eady problem (Eady, 1949). In the Eady problem, perturbations evolve around a basic state with constant shear $\Lambda$ in thermal wind balance with a constant meridional buoyancy gradient $f\Lambda$, such that

$$u = \underbrace{\Lambda z}_{\stackrel{\text{def}}{=} U} + u' \,, \quad \text{and} \quad b = \underbrace{-f\Lambda y}_{\stackrel{\text{def}}{=} B} + b' \,. \tag{10}$$

We use Oceananigans' BackgroundFields to simulate the nonlinear evolution of $(u', v, w)$ and $b'$ expanded around $U$ and $B$ in a doubly-periodic domain. We impose an initially stable density stratification with $b' = N^2 z$ and $N^2 = 10^{-7}$ s$^{-2}$ superposed with random noise. The Richardson number of the initial condition is $Ri = N^2/\partial_z U = N^2/\Lambda$; we choose mean shear $\Lambda$ so that $Ri = 1$, which guarantees the basic state is unstable to baroclinic instability but stable to symmetric and Kelvin-Helmholtz instability (Stone, 1971). A portion of the script is shown in listing 6.

**Listing 6:** Large eddy simulation of the Eady problem expanded around the background geostrophic shear with Ri = 1.

```julia
grid = RectilinearGrid(GPU(); size = (1024, 1024, 64),
    x = (0, 4096), y = (0, 4096), z = (0, 128),
    topology = (Periodic, Periodic, Bounded))

f, N2, Ri = 1e-4, 1e-7, 1
parameters = (f=f, Lambda=sqrt(N2/Ri)) # U = Lambda*z, so Ri = N2 / dz(U) = N2 / Lambda

@inline U(x, y, z, t, p) = + p.Lambda * z
@inline B(x, y, z, t, p) = - p.f * p.Lambda * y

background_fields = (u = BackgroundField(U; parameters),
    b = BackgroundField(B; parameters))

model = NonhydrostaticModel(; grid, background_fields,
    advection = WENO(order=9), coriolis = FPlane(; f),
    tracers = :b, buoyancy = BuoyancyTracer())

dz = minimum_zspacing(grid)
bi(x, y, z) = N2 * z + 1e-2 * N2 * dz * (2rand() - 1)
set!(model, b=bi)
```

Our Eady simulation uses fully-turbulence-resolving resolution with 4 meter horizontal spacing and 2 meter vertical spacing in a 4 km x 4 km x 128 m domain and simulates 30 days on a single Nvidia H100 GPU. Four snapshots of vertical vorticity normalized by $f$ (the Rossby number) are shown in figure 4, illustrating the growth of kilometer-scale vortex motions amid bursts of meter-scale three-dimensional turbulence that develop along thin filaments of vertical vorticity and vertical shear. This configuration captures a competition between baroclinic instability, which acts to "restratify" or strengthen boundary layer stratification, and three-dimensional turbulent mixing driven by a forward cascade from kilometer-scale motions (Molemaker et al., 2010; Dong et al., 2024). Shear and convection associated with atmospheric storms provide an additional important source of ocean turbulent mixing interacting with kilometer-scale baroclinic instability (Boccaletti et al., 2007; Callies & Ferrari, 2018); these additional processes can be simulated by adding surface fluxes to our Eady configuration. Future work is required to understand what resolution is required to faithfully simulate multi-scale flows such as Eady turbulence, in the presence of realistically-strong density stratification and when using high-order WENO advection schemes.

[Figure 4: Surface vertical vorticity in a large eddy simulation of the Eady problem with Ri = 1 initially, after t = 4.6, 6, 7.7, and 20 days. The grid spacing is 4 x 4 x 2 meters in x, y, z. Part of the script that produces this simulation is shown in listing 6.]

Finally, we illustrate Oceananigans' capabilities for realistic, three-dimensional large eddy simulations in complex geometries by simulating temperature- and salinity-stratified tidal flow past a headland, reminiscent of an extensively observed and modeled flow past Three Tree Point in Puget Sound in the Pacific Northwest of the United States (Pawlak et al., 2003; Warner & MacCready, 2014). The bathymetry involves a sloping wedge that juts from a square-sided channel, such that

$$z_b(x, y) = -H\left(1 + \frac{y + |x|}{\delta}\right) \,, \tag{11}$$

where $\delta = L/2$ represents the scale of the bathymetry, $L$ is the half-channel width in $y$ (the total width is $2L$), and $H = 128$ m is the depth of the channel, and $z = z_b(x, y)$ is the height of the bottom. The flow is driven by a tidally-oscillating boundary velocity

$$U(t) = U_2 \sin\left(\frac{2\pi t}{T_2}\right) \tag{12}$$

imposed at the east and west boundaries. Here, $T_2 = 12.421$ hours is the period of the semi-diurnal lunar tide, and $U_2 = 0.15$ m s$^{-1}$ is the characteristic tidal velocity around Three Tree Point. The initial temperature and salinity are

$$T|_{t=0} = 12 + 4\frac{z}{H} \; ^\circ\text{C} \,, \quad \text{and} \quad S|_{t=0} = 32 \; \text{g kg}^{-1} \,. \tag{13}$$

A portion of the script that implements this simulation is shown in listing 7.

**Listing 7:** Large eddy simulation of flow past a headland reminiscent of Three Tree Point in the Pacific Northwest (see Pawlak et al., 2003; Warner & MacCready, 2014).

```julia
H, L = 256meters, 1024meters
delta = L / 2
x, y, z = (-3L, 3L), (-L, L), (-H, 0)
Nz = 64

grid = RectilinearGrid(GPU(); size=(6Nz, 2Nz, Nz), halo=(6, 6, 6),
    x, y, z, topology=(Bounded, Bounded, Bounded))

wedge(x, y) = -H * (1 + (y + abs(x)) / delta)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(wedge))

T2 = 12.421hours
U2 = 0.1 # m/s

@inline U(x, y, z, t, p) = p.U2 * sin(2pi * t / p.T2)
@inline U(y, z, t, p) = U(zero(y), y, z, t, p)

open_bc = OpenBoundaryCondition(U; parameters=(; U2, T2),
    scheme = PerturbationAdvection(inflow_timescale = 2minutes,
                                   outflow_timescale = 2minutes))

u_bcs = FieldBoundaryConditions(east = open_bc, west = open_bc)

@inline ambient_temperature(x, z, t, H) = 12 + 4z/H
@inline ambient_temperature(x, y, z, t, H) = ambient_temperature(x, z, t, H)
ambient_temperature_bc = ValueBoundaryCondition(ambient_temperature; parameters = H)
T_bcs = FieldBoundaryConditions(east = ambient_temperature_bc,
                                west = ambient_temperature_bc)

ambient_salinity_bc = ValueBoundaryCondition(32)
S_bcs = FieldBoundaryConditions(east = ambient_salinity_bc, west = ambient_salinity_bc)

model = NonhydrostaticModel(; grid, tracers = (:T, :S),
    buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState()),
    advection = WENO(order=9), coriolis = FPlane(latitude=47.5),
    boundary_conditions = (; T=T_bcs, u = u_bcs, S = S_bcs))

Ti(x, y, z) = ambient_temperature(x, y, z, 0, H)
set!(model, T=Ti, S=32)
```

The oscillatory, turbulent flow is visualized in figure 5. The calculation of Ertel Potential Vorticity shown in figure 5c uses the companion package Oceanostics (Chor et al., 2025). This simulation uses a grid fitted immersed boundary method, which produces a "staircase" representation of the bathymetry. Future work remains to develop a cut cell method (e.g. Adcroft et al., 1997; Yamazaki et al., 2016) to produce a smoother, piecewise linear representation of bathymetry, which will prevent the generation of noise and spurious waves occurring at the sharp corners of staircase bathymetry.

[Figure 5: Along-channel velocity, temperature, and Ertel potential vorticity in a tidally-forced flow past an idealized headland with open boundaries. The tidal flow occurs in the x-directions and the snapshot depicts the flow just after the tide has turned to the negative-x direction.]

### 3.2 Hydrostatic model with a free surface

The HydrostaticFreeSurfaceModel solves the hydrostatic, rotating Boussinesq equations with a free surface. The hydrostatic approximation, inherent to the HydrostaticFreeSurfaceModel, means that the vertical momentum equation used by NonhydrostaticModel, $\hat{z} \cdot$ (6), is replaced by a statement of hydrostatic balance,

$$\partial_z p = b \,, \tag{14}$$

while the vertical velocity is obtained diagnostically from the continuity equation,

$$\partial_z w = -\nabla_h \cdot \mathbf{u}_h \,. \tag{15}$$

As a result, time-stepping the HydrostaticFreeSurfaceModel does not require solving a three-dimensional Poisson equation for pressure. Moreover, the HydrostaticFreeSurfaceModel introduces a free surface displacement $\eta$, which obeys the linearized equation

$$\partial_t \eta = w|_{z=0} \,. \tag{16}$$

Equation (16) replaces the rigid-lid impenetrability condition $w|_{z=0} = 0$ typically applied at top boundaries in the NonhydrostaticModel. The numerical algorithms and computational performance of the HydrostaticFreeSurfaceModel are described in more detail by Silvestri et al. (2025).

In the HydrostaticFreeSurfaceModel, the horizontal momentum $\mathbf{u}_h = u\,\hat{x} + v\,\hat{y}$ evolves according to

<!-- [Equation 17: HydrostaticFreeSurfaceModel horizontal momentum equation] -->
$$\partial_t \mathbf{u}_h = -\nabla_h p \underbrace{- g\nabla_h \eta}_{\text{free surface}} \underbrace{- (\mathbf{u} \cdot \nabla)\mathbf{u}_h}_{\text{momentum advection}} \underbrace{- \mathbf{f} \times \mathbf{u}}_{\text{Coriolis}} \underbrace{- \nabla \cdot \boldsymbol{\tau}}_{\text{closure}} + \underbrace{\mathbf{F}_{u_h}}_{\text{forcing}} \,, \tag{17}$$

where $p$ is the hydrostatic kinematic pressure anomaly, $\eta$ is the free surface displacement, $\mathbf{u} = u\,\hat{x} + v\,\hat{y} + w\,\hat{z}$ is the three-dimensional velocity, $\mathbf{f}$ is the background vorticity associated with a rotating frame of reference, $\boldsymbol{\tau}$ is the stress associated with subgrid turbulent horizontal momentum transport, and $\mathbf{F}_{u_h}$ is a body force. As described in appendix A, equation (17) may be integrated with a predictor-corrector split-explicit method that uses short substeps for the "barotropic mode" (the fast component of the free surface displacement $\eta$ and vertically-integrated horizontal momentum) and a much slower step for the three-dimensional part of the solution. Implicit and fully explicit methods for time-integration of (A16)--(17) are also available.

Momentum advection in HydrostaticFreeSurfaceModel can be formulated in three ways: in the same flux form used by NonhydrostaticModel,

$$(\mathbf{u} \cdot \nabla)\mathbf{u}_h = \nabla \cdot (\mathbf{u}\,\mathbf{u}_h) \,, \tag{18}$$

in a standard "vector invariant" form that facilitates simulations on curvilinear grids (Adcroft et al., 2004),

$$(\mathbf{u} \cdot \nabla)\mathbf{u}_h = \zeta\,\hat{z} \times \mathbf{u}_h + w\,\partial_z \mathbf{u}_h + \nabla_h \tfrac{1}{2}|\mathbf{u}_h|^2 \,, \tag{19}$$

and a modified vector invariant form,

$$(\mathbf{u} \cdot \nabla)\mathbf{u}_h = \zeta\,\hat{z} \times \mathbf{u}_h - \mathbf{u}_h\,\partial_z w + \partial_z(w\,\mathbf{u}_h) + \nabla_h \tfrac{1}{2}|\mathbf{u}_h|^2 \,, \tag{20}$$

that is used for a WENO-based vector invariant scheme (Silvestri et al., 2024). The WENO-based vector invariant scheme is particularly well-suited for simulating geophysical turbulence subject to strong Coriolis forces, because selectively dissipates enstrophy and the variance of divergence and does not require an explicit turbulence closure for numerical stability (Silvestri et al., 2024). (Zhang et al. (2025) show how the WENO vector invariant methodology may be translated to apply to potential vorticity in a layered formulation, and thereby dissipate potential vorticity.) This contrasts the WENO vector invariant scheme with second-order flux form or standard vector invariant schemes, which produce oscillatory errors and must be paired with an explicit turbulence closure, and with a flux form WENO scheme, which targets the dissipation of kinetic energy rather than enstrophy.

Tracer evolution in the HydrostaticFreeSurfaceModel is governed by the conservation law

<!-- [Equation 21: HydrostaticFreeSurfaceModel tracer conservation] -->
$$\partial_t c = \underbrace{-\nabla \cdot [(\mathbf{u} + \mathbf{u}_p)\,c] + c\,\nabla \cdot \mathbf{u}_p}_{\text{tracer advection}} \underbrace{- \nabla \cdot \mathbf{J}_c}_{\text{closure}} + \underbrace{S_c}_{\text{biogeochemistry}} + \underbrace{F_c}_{\text{forcing}} \,, \tag{21}$$

which is identical to the NonhydrostaticModel tracer equation (9), except that background fields are not supported. HydrostaticFreeSurfaceModel also supports prescribing the velocity field $\mathbf{u}$ in (21), enabling cheap tracer advection and Lagrangian particle trajectory simulations.

**Listing 8:** Two-dimensional simulation of tidally-forced stratified flow over a superposition of randomly-positioned Gaussian seamounts. Results are shown in Figure 6.

```julia
using Oceananigans, Oceananigans.Units

grid = RectilinearGrid(size = (2000, 200), halo = (4, 4),
    x = (-1000kilometers, 1000kilometers),
    z = (-2kilometers, 0),
    topology = (Periodic, Flat, Bounded))

h0 = 100          # typical mountain height (m)
delta = 20kilometers # mountain width (m)
seamounts = 42
W = grid.Lx - 4delta
x0 = W .* (rand(seamounts) .- 1/2) # mountains' positions in [-Lx/2+2delta, Lx/2-2delta]
h  = h0 .* (1 .+ rand(seamounts))  # mountains' heights in [h0, 2h0]

bottom(x) = -grid.Lz + sum(h[s] * exp(-(x - x0[s])^2 / 2delta^2) for s = 1:seamounts)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom))

T2 = 12.421hours # period of M2 tide constituent
@inline tidal_forcing(x, z, t, p) = p.U2 * 2pi / p.T2 * sin(2pi / p.T2 * t)
u_forcing = Forcing(tidal_forcing, parameters=(; U2=0.1, T2=T2))

model = HydrostaticFreeSurfaceModel(; grid, tracers=:b, buoyancy=BuoyancyTracer(),
    momentum_advection=WENO(), tracer_advection=WENO(),
    forcing=(; u=u_forcing))

bi(x, z) = 1e-5 * z
set!(model, b=bi)
```

Listing 8 implements a two-dimensional simulation of tidally-forced stratified flow over a series of randomly-positioned Gaussian seamounts, using flux form WENO advection schemes to dissipate kinetic energy and tracer variance, as appropriate for a two-dimensional simulation dominated by wave motion. In this simulation, the horizontal flow associated with tidal forcing is forced upwards over the seamounts, perturbing the stable, vertically-stratified buoyancy distribution and forcing the development of internal inertia-gravity waves, sometimes called "internal tides" (Garrett & Kunze, 2007). The vertical velocity is visualized in figure 6, showing that the tidal signal is dominated by "mode-one", whose structure extends smoothly from the top to the bottom of the domain.

[Figure 6: Vertical velocity of an internal wave field excited by tidally-forced stratified flow over a superposition of randomly-positioned Gaussian seamounts, after 16 tidal periods.]

#### 3.2.1 Vertical mixing and mesoscale turbulence parameterizations

Similar to the NonhydrostaticModel, the closure tracer fluxes $\mathbf{J}_c$ may be three-dimensionally, horizontally, or vertically Laplacian or biharmonic. In addition, a closure is provided that computes both isopycnal "symmetric" diffusive tracer fluxes (Redi, 1982) and skew diffusive tracer fluxes (Gent & Mcwilliams, 1990; Griffies, 1998), which are typically used to parameterize stirring and vertical form stress, respectively, associated with unresolved mesoscale turbulence (e.g. Young, 2012). The symmetric contribution (often called "Redi") and skew contribution (often called "GM" after Gent--McWilliams) to the tracer flux may be tapered in the presence of a large isopycnal slope. The GM flux may be formulated as tracer advection, permitting combination with the other velocity components $\mathbf{u}$ and $\mathbf{u}_p$, but with the downside of invoking additional gradients of the isopycnal slope (Griffies, 1998).

Oceananigans provides a number of closures that parameterize vertical fluxes of tracers and momentum associated with unresolved turbulence with scales from $O(1)$--$O(100)$ meters. Depending on the parameterization, the evolution of auxiliary tracers like turbulent kinetic energy and the turbulent kinetic energy dissipation rate may also be simulated. Vertical mixing parameterizations are useful for hydrostatic simulations where vertical mixing is otherwise unresolved due to a coarse horizontal grid spacing. For example, such regional and global configurations, horizontal grid spacing typically varies from $O(100\;\text{m})$ to $O(100\;\text{km})$.

**Listing 9:** Comparison of two vertical mixing parameterizations in the evolution of an initially linearly stratified boundary layer subjected to stationary surface fluxes of buoyancy and momentum. Results are shown in Figure 7.

```julia
using Oceananigans
using Oceananigans.Units

function vertical_mixing_simulation(closure; N2=1e-5, Jb=1e-7, tx=-5e-4)
    grid = RectilinearGrid(size=50, z=(-200, 0), topology=(Flat, Flat, Bounded))
    buoyancy = BuoyancyTracer()

    b_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(Jb))
    u_bcs = FieldBoundaryConditions(top=FluxBoundaryCondition(tx))

    model = HydrostaticFreeSurfaceModel(; grid, closure, tracers=:b, buoyancy,
        boundary_conditions=(u=u_bcs, b=b_bcs))

    bi(z) = N2 * z
    set!(model, b=bi)

    simulation = Simulation(model, dt=1minute, stop_time=24hours)
    return run!(simulation)
end
```

Listing 9 implements a simulation of wind-driven vertical mixing in a single column model using two parameterizations: CATKE (Wagner, Hillier, et al., 2025), which has one additional equation for the evolution of turbulent kinetic energy (TKE), and k-epsilon (Umlauf & Burchard, 2005), which has two additional equations for both TKE and TKE dissipation. Figure 7 plots the result, showing how k-epsilon mixes less than CATKE. This discrepancy in mixing rates is likely due to differences in how the models are calibrated. While all of CATKE's parameters are jointly calibrated to 35 large eddy simulations (LES) that include surface wave effects (Wagner, Hillier, et al., 2025), k-epsilon parameters are calibrated one-by-one by referencing laboratory experiments and observations of increasing complexity (Umlauf & Burchard, 2003). An Ri-based scheme similar to the one proposed by Pacanowski and Philander (1981) is also available. Both CATKE and k-epsilon may be optionally substepped, a useful performance optimization for coarse resolution simulations that can otherwise accommodate relatively long baroclinic time steps. Directions for future work in parameterization include the implementation of established closures for Langmuir turbulence (for example, Reichl & Li, 2019; Harcourt, 2015), the development of new closures (for example following Legay et al., 2025; Wagner & Constantinou, 2025), calibration of closures like k-epsilon following the approach by Wagner, Hillier, et al. (2025), and the implementation of energy-constrained formulations of the GM and Redi tracer flux (Mak et al., 2018; Jansen et al., 2019).

[Figure 7: Results from two vertical mixing parameterizations: CATKE and k-epsilon, implemented as described in Listing 9.]

#### 3.2.2 Ocean simulations on global and near-global spherical shell grids

The HydrostaticFreeSurfaceModel can be used to simulate regional or global ocean circulation on rectilinear grids, latitude-longitude grids, and the tripolar grid (Murray, 1996) to cover the entirety of Earth's global ocean with the current continental configuration. A cubed sphere grid that covers the entire sphere and is therefore useful for aquaplanet simulations, is also implemented and currently being validated. To illustrate simulations on global and near-global grids, we simulate baroclinic instability on three spherical grids: a latitude-longitude grid, a tripolar grid with islands placed over the north pole singularities, and a "rotated" latitude-longitude grid wherein the polar grid singularities are rotated from the geographic north pole to 55 N, 70 W, corresponding to the default location of the tripolar grid's two north poles. A portion of the code that produces the simulation is given in listing 10 and the results are visualized in figure 8. In addition to completing the validation of cubed sphere simulations, an interesting direction for future work is to develop capabilities for aquaplanet simulations on latitude-longitude grids, using similar approaches as used for atmosphere models (Richardson et al., 2007; Khairoutdinov et al., 2022). All of Oceananigans' grids currently use a shallow or "thin shell" approximation; relaxing this approximation is another direction for future work.

**Listing 10:** Near-global simulations of baroclinic instability on three different grids.

```julia
arch = GPU()
Nx, Ny, Nz = size = (4 * 360, 4 * 170, 10)
halo = (7, 7, 7)
H = 3000
latitude, longitude, z = (-85, 85), (0, 360), (-H, 0)

lat_lon_grid = LatitudeLongitudeGrid(arch; size, halo, latitude, longitude, z)
rotated_grid = RotatedLatitudeLongitudeGrid(arch; size, halo, latitude,
    longitude, z, north_pole=(70, 55))

# TripolarGrid with "Gaussian islands" over the two north poles
underlying_tripolar_grid = TripolarGrid(arch; size, halo, z)

dphi, dlambda, lambda0, phi0 = 4, 8, 70, 55
isle(lambda, phi) = ((lambda - lambda0)^2 / 2dlambda^2 + (phi - phi0)^2 / 2dphi^2) < 1
cylindrical_isles(lambda, phi) = H * (isle(lambda, phi) + isle(lambda - 180, phi) - 1)
tripolar_grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(cylindrical_isles))

momentum_advection = WENOVectorInvariant(order=9)
tracer_advection = WENO(order=7)
coriolis = HydrostaticSphericalCoriolis()
buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())
free_surface = SplitExplicitFreeSurface(grid, substeps=60)

grid = lat_lon_grid # rotated_grid, tripolar_grid
model = HydrostaticFreeSurfaceModel(; grid, momentum_advection, tracer_advection,
    coriolis, free_surface, buoyancy, tracers=(:T, :S))

Ti(lambda, phi, z) = 30 * (1 - tanh((abs(phi) - 45) / 8)) / 2 + rand()
Si(lambda, phi, z) = 28 - 5e-3 * z + rand()
set!(model, T=Ti, S=Si)

simulation = Simulation(model, dt=2minutes, stop_time=180days)
```

[Figure 8: Visualization of relative vorticity (top row) and temperature (bottom row) in simulations of baroclinic instability on three different grids: a latitude-longitude grid (left), a tripolar grid with cylindrical islands centered on the two north pole singularities (middle), and a "rotated" latitude-longitude grid with a grid pole located at 55 N, 70 W. Listing 10.]

#### 3.2.3 Realistic ocean simulations beneath prescribed atmospheric states with ClimaOcean

The coupled modeling package ClimaOcean (Wagner, Silvestri, et al., 2025) implements a framework for coupled modeling including model integration, diagnostics, flux, interfacial state computations. In ClimaOcean, turbulent interfacial fluxes are computed using Monin--Obukhov similarity theory (Monin & Obukhov, 1954) following Edson et al. (2014) for air-sea fluxes and Grachev et al. (2007) for air-ice fluxes. ClimaOcean additionally provides utilities for downloading and interfacing with JRA55 reanalysis data (Tsujino et al., 2018), building grids based on Earth bathymetry and initializing simulations from datasets including the Estimating the Circulation and Climate of the Ocean (ECCO) state estimate (Forget et al., 2015), Mercator Ocean's Global Ocean Reanalysis (GLORYS), and EN4 (Good et al., 2013).

We illustrate ClimaOcean's capabilities by implementing a coupled ocean and sea ice simulation on a tripolar grid with 1/6th degree resolution. The ocean component is based on Oceananigans' HydrostaticFreeSurfaceModel, using 9th-order WENOVectorInvariant advection for momentum and 7th order WENO advection for tracers, CATKE for vertical mixing, the TEOS10EquationOfState, and bathymetry derived from ETOPO1. The sea ice component is based on ClimaSeaIce, using a 0-layer "slab" thermodynamic formulation with constant heat conductivity and a viscoplastic rheology implemented with the Kimmritz et al. (2016) solver. Part of our code, which uses a prescribed atmospheric state and radiation fields derived from the JRA55 reanalysis, and initializes ocean temperature, salinity, sea ice thickness, and sea ice concentration from the ECCO state estimate, is shown in listing 11.

**Listing 11:** A coupled ocean and sea-ice simulation beneath a prescribed atmospheric state from the JRA55 reanalysis (Tsujino et al., 2018) using ClimaOcean. Ocean temperature, ocean salinity, sea ice thickness, and sea ice concentration are initialized using the ECCO state estimate (Forget et al., 2015).

```julia
Nx, Ny, Nz = 2160, 1080, 60 # 1/6th degree
z = ExponentialDiscretization(Nz, -6000, 0, mutable=true) # enables z* vertical coordinate
underlying_grid = TripolarGrid(GPU(); z, size=(Nx, Ny, Nz), halo=(7, 7, 7))
bathymetry = ClimaOcean.regrid_bathymetry(grid) # based on ETOPO1
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bathymetry))

# Build an ocean and sea ice simulations and initialize
# to the ECCO state estimate on Jan 1, 1993
ocean = ClimaOcean.ocean_simulation(grid)
sea_ice = ClimaOcean.sea_ice_simulation(grid)

date = CFTime.DateTimeProlepticGregorian(1993, 1, 1)
set!(ocean.model, T = ClimaOcean.ECCOMetadata(:temperature; date),
                  S = ClimaOcean.ECCOMetadata(:salinity; date))

set!(sea_ice.model, h = ClimaOcean.ECCOMetadata(:sea_ice_thickness; date),
                    a = ClimaOcean.ECCOMetadata(:sea_ice_concentration; date))

# Forced by JRA55 reanalysis
backend = ClimaOcean.JRA55NetCDFBackend(41)
atmosphere = ClimaOcean.JRA55_prescribed_atmosphere(arch; backend)

coupled_model = ClimaOcean.OceanSeaIceModel(ocean, sea_ice; atmosphere)
```

The ocean surface speed (blue to yellow colors) and sea ice speed (red to white colors) generated after 180 days of simulation time is shown in figure 9. For more information about Oceananigans performance on multiple GPUs, see Silvestri et al. (2025).

[Figure 9: Ocean surface and ice speed (masked for ice concentrations less than 10^-3) in a coupled ocean and sea-ice simulation on a 1/6th degree tripolar grid (Murray, 1996) from two perspectives in space. The ocean and sea ice state is initialized from the ECCO state estimate (Forget et al., 2015) and forced by the JRA55-do atmospheric reanalysis (Tsujino et al., 2018). Air-sea fluxes and air-ice fluxes and skin temperature are computed by fixed point iteration of bulk formula. The air-sea bulk formulation uses Edson et al. (2014) stability functions. The air-ice bulk formulation uses the Grachev et al. (2007) and Paulson (1970) stability functions for stable and unstable conditions, respectively. The sea ice model uses a 0-layer thermodynamic formulation with constant heat conductivity and a viscoplastic rheology with the solver by Kimmritz et al. (2016). The simulation is spun up with a 20 second time-step for 60 days and then integrated for 58 years with a 6 minute time step.]

## 4 Conclusions

This paper describes the GPU-based ocean modeling software called "Oceananigans" written in the high-level Julia programming language. Oceananigans provides a productive script-based user interface and reduces the cost of high-resolution simulations of oceanic motion at any scale. The current state of Oceananigans realizes a particular strategy for improving dynamical cores: basic C-grid WENO numerics for turbulence resolving simulations coupled to the raw power of GPU acceleration.

Because Oceananigans enables high-resolution simulations with resources as small as a single GPU, it increases access to high-fidelity ocean modeling. But it also enables a new class of ultra-high-resolution simulations. For example, on the Perlmutter supercomputer (National Energy Research Scientific Computing Center, 2025), it is currently possible to complete a 100-member ensemble of century-long global ocean simulations at 10 kilometer resolution in 10 days of wall time -- thereby resolving mesoscale turbulent mixing, a prominent bias in ocean models and a fundamental process missing from most climate simulations today. These new capabilities address uncertainty in ocean heat and carbon uptake in climate projections.

Oceananigans user interface facilitates composition with other software packages and has fostered the development of an ecosystem of packages for ocean modeling. For example, OceanBioME (Strong-Wright et al., 2023) implements Oceananigans-compatible biogeochemistry models, oriented towards ecosystem dynamics and compatible with both the hydrostatic and nonhydrostatic models. A second biogeochemistry package is also under development for climate simulations. The Oceanostics (Chor et al., 2025) package implements complex diagnostics in Oceananigans syntax, useful for online and offline analysis of large eddy simulations.

Initial work is also under way to couple Oceananigans-based ocean models with prognostic atmosphere models, including the Climate Modeling Alliance atmosphere dynamical core (Yatunin et al., 2025) and the simpler SpeedyWeather (Klower et al., 2024). A second effort called DJ4Earth is using Enzyme (Moses et al., 2021) and Reactant to develop an adjoint for Oceananigans, and to more generally use auto-differentiation to compute the gradients of cost functions that invoke Oceananigans simulations. Additional directions for future work include:

- Implementing, developing, and calibrating energy-based parameterizations for mesoscale turbulence (e.g. Mak et al., 2018; Jansen et al., 2019);
- Implementing, developing, and calibrating more vertical mixing parameterizations including for Langmuir turbulence (for example, building upon Reichl & Li, 2019; Harcourt, 2015; Legay et al., 2025; Wagner & Constantinou, 2025);
- Characterization of numerical mixing and dissipation in simulations with WENO-based advection schemes, especially for long running climate simulations;
- Development of ClimaSeaIce, including a multi-layer, multi-category thermodynamic component and elasto-brittle rheology (Dansereau et al., 2016);
- The development and calibration of parameterizations for the air-sea interface state, air-ice interface state, and coupling (e.g. Pelletier et al., 2021);
- Further performance optimization especially in multi-GPU configurations, using bespoke GPU features and mixed precision to further optimize GPU performance;
- Development of new numerical methods for representing bathymetry (e.g. shaved cells, Adcroft et al., 1997; Shaw & Weller, 2016; Yamazaki et al., 2016);
- Generation of high-resolution coupled ocean and sea ice training data for training AI ocean emulators (e.g. Dheeshjith et al., 2025);
- Development and coupling with a GPU-native surface model and coupling with PiCLES, a recently-developed Particle-In-Cell surface wave model (Hell et al., 2024);
- Development of global coupled atmosphere-ocean configurations with biogeochemistry;
- Development of open boundary conditions for HydrostaticFreeSurfaceModel for regional configurations similar to the implementation used for NonhydrostaticModel;
- Coupling with hybrid physics/ML atmosphere models (Kochkov et al., 2024) and fully-ML atmosphere models like ACE (Watt-Meyer et al., 2023, 2024) and GraphCast (Lam et al., 2023).

Driven by the recent explosion in computationally-intensive applications of machine learning (ML), computational science and technology is advancing at an accelerating rate. To keep pace -- to continue to use the world's fastest computers, to maintain scientific productivity commensurate with other fields, and to enable the next generation of theory and ML-based parameterizations -- the development of ocean modeling software must also accelerate. Oceananigans represents a step towards the faster development of ocean modeling software. But substantial and sustained growth in the developer community is still required for Oceananigans to reach its potential.

## Appendix A: Time stepping and time discretization

In this section we describe time stepping methods and time discretization options for the NonhydrostaticModel and the HydrostaticFreeSurfaceModel.

### A1 Time discretization for tracers

Tracers are stepped forward with similar schemes in the NonhydrostaticModel and the HydrostaticFreeSurfaceModel, each of which includes optional implicit treatment of vertical diffusion terms. Equation (9) is abstracted into two components,

$$\partial_t c = G_c + \partial_z(\kappa_z \partial_z c) \,, \tag{A1}$$

where, if specified, $\kappa_z$ is the vertical diffusivity of $c$ to be treated with a VerticallyImplicitTimeDiscretization, and $G_c$ is the remaining component of the tracer tendency from equation 9. (Vertical diffusion treated with an ExplicitTimeDiscretization is also absorbed into $G_c$.) We apply a semi-implicit time discretization of vertical diffusion to approximate integral of (A1) from $t^m$ to $t^{m+1}$,

$$(1 - \Delta t\, \partial_z \kappa_z^m \partial_z)\, c^{m+1} = c^m + \int_{t^m}^{t^{m+1}} G_c \, dt \,, \tag{A2}$$

where $\Delta t \stackrel{\text{def}}{=} t^{m+1} - t^m$. The tendency integral $\int_{t^m}^{t^{m+1}} G_c \, dt$ is evaluated either using a "quasi"-second order Adams-Bashforth scheme (QAB2, which is formally first-order), or a low-storage third-order Runge-Kutta scheme (RK3). For QAB2, the integral in (A2) spans the entire time-step and takes the form

$$\frac{1}{\Delta t} \int_{t^m}^{t^{m+1}} G_c \, dt \approx \left(\frac{3}{2} + \chi\right) G_c^m - \left(\frac{1}{2} + \chi\right) G_c^{m-1} \,, \tag{A3}$$

where $\chi$ is a small parameter, chosen by default to be $\chi = 0.1$. QAB2 requires one tendency evaluation per time-step. For RK3, the indices $m = (1, 2, 3)$ correspond to substages, and the integral in (A2) takes the form

$$\frac{1}{\Delta t} \int_{t^m}^{t^{m+1}} G_c \, dt \approx \gamma_m G_c^m - \zeta_m G_c^{m-1} \,, \tag{A4}$$

where $\gamma = (8/15, 5/12, 3/4)$ and $\zeta = (0, 17/60, 5/12)$ for $m = (1, 2, 3)$ respectively. RK3 requires three evaluations of the tendency $G_c$ per time-step. RK3 is self-starting because $\zeta_1 = 0$, while QAB2 must be started with a forward-backwards Euler step (the choice $\chi = -1/2$ in (A3)). Equation (A2) is solved with a tridiagonal algorithm following a second-order spatial discretization of $\partial_z \kappa_z^n \partial_z c^{m+1}$ -- either once per time-step for QAB2, or three times for each of the RK3's three stages.

VerticallyImplicitTimeDiscretization permits longer time-steps when using fine vertical spacing. Listing 12 illustrates the differences between vertically-implicit and explicit time discretization using one-dimensional diffusion by a top-hat diffusivity profile. The results are shown in figure A1.

**Listing 12:** Diffusion of a tracer by a top hat tracer diffusivity profile using various time steps and time discretizations.

```julia
using Oceananigans

grid = RectilinearGrid(size=20, z=(-2, 2), topology=(Flat, Flat, Bounded))
time_discretization = VerticallyImplicitTimeDiscretization()
kappa(z, t) = exp(-z^2)
closure = VerticalScalarDiffusivity(time_discretization; kappa)
model = HydrostaticFreeSurfaceModel(; grid, closure, tracers=:c)
```

[Figure A1: Simulations of tracer diffusion by a top hat diffusivity profile using various choices of time-discretization and time-step size. With a long time-step of dt = 0.5, ExplicitTimeDiscretization is unstable while VerticallyImplicitTimeDiscretization is stable. Let the vertically-implicit solution depends on the long time-step dt = 0.5, as revealed by comparison with ExplicitTimeDiscretization using dt = 10^-4.]

### A2 The pressure correction method for momentum in NonhydrostaticModel

The NonhydrostaticModel uses a pressure correction method for the momentum equation (6) that ensures $\nabla \cdot \mathbf{u} = 0$. We rewrite (6) as

$$\partial_t \mathbf{u} = -\nabla p + b\,\hat{z} + \mathbf{G}_u + \partial_z(\nu_z \partial_z \mathbf{u}) \,, \tag{A5}$$

where, if specified, $\nu_z$ is the vertical component of the viscosity that will be treated with a vertically-implicit time discretization, $\nabla p$ is the total pressure gradient, and $\mathbf{G}_u$ is the rest of the momentum tendency. We decompose $p$ into a "hydrostatic anomaly" $p'$ tied to the density anomaly $\rho'$, and a nonhydrostatic component $p_{nh}$, such that

$$p = p_{nh} + p_{hy} \,, \quad \text{where} \quad \partial_z p_{hy} \stackrel{\text{def}}{=} b \,. \tag{A6}$$

By computing $p_{hy}$ in (A6), we recast (A5) without $b$ and with $\nabla p = \nabla p_{nh} + \nabla_h p_{hy}$. Next, integrating (A5) in time from $t^m$ to $t^{m+1}$ yields

$$\mathbf{u}^{m+1} = \mathbf{u}^m + \int_{t^m}^{t^{m+1}} [\mathbf{G}_u - \nabla p_{nh} + \partial_z(\nu_z \partial_z \mathbf{u})] \, dt \,. \tag{A7}$$

Next we introduce the predictor velocity $\tilde{\mathbf{u}}$, defined such that

$$(1 - \Delta t\, \partial_z \nu_z^m \partial_z)\, \tilde{\mathbf{u}} = \mathbf{u}^m + \int_{t^m}^{t^{m+1}} \mathbf{G}_u \, dt \,, \tag{A8}$$

or in other words, defined as a velocity-like field that cannot feel nonhydrostatic pressure gradient $\nabla p_{nh}$. Equation (A8) uses a semi-implicit treatment of vertical momentum diffusion which is similar but slightly different to the treatment of tracer diffusion in (A2),

$$\int_{t^m}^{t^{m+1}} \partial_z(\nu_z \partial_z \mathbf{u}) \, dt \approx \Delta t\, \partial_z(\nu_z^m \partial_z \tilde{\mathbf{u}}) \,. \tag{A9}$$

The integral in (A8) is evaluated with the same methods used for tracers -- either (A3) for QAB2 or (A4) when using RK3. With a second-order discretization of vertical momentum diffusion, the predictor velocity in (A8) may be computed with a tridiagonal solver.

Introducing a fully-implicit time discretization for $p_{nh}$,

$$\int_{t^m}^{t^{m+1}} \nabla p_{nh} \, dt \approx \Delta t\, \nabla p_{nh}^{m+1} \,, \tag{A10}$$

and inserting (A10) into (A8), we derive the pressure correction to the predictor velocity,

$$\mathbf{u}^{m+1} - \tilde{\mathbf{u}} = -\Delta t\, \nabla p_{nh}^{m+1} \,. \tag{A11}$$

The final ingredient needed to complete the pressure correction scheme is an equation for the nonhydrostatic pressure $p_{nh}^{m+1}$. For this we form $\nabla \cdot$ (A11) and use $\nabla \cdot \mathbf{u}^{m+1} = 0$ to obtain a Poisson equation for $p_{nh}^{m+1}$,

$$\nabla^2 p_{nh}^{m+1} = \frac{\nabla \cdot \tilde{\mathbf{u}}}{\Delta t} \,. \tag{A12}$$

Boundary conditions for equation (A12) may be derived by evaluating $\hat{n} \cdot$ (A7) on the boundary of the domain.

On RectilinearGrids, we solve (A12) using an eigenfunction expansion of the discrete second-order Poisson operator $\nabla^2$ evaluated via the fast Fourier transform (FFT) in equispaced directions (Schumann & Sweet, 1988) plus a tridiagonal solve in variably-spaced directions. With the FFT-based solver, boundary conditions on $p_{nh}^{m+1}$ are accounted for by enforcing $\hat{n} \cdot \tilde{\mathbf{u}} = \hat{n} \cdot \mathbf{u}^{m+1}$ on boundary cells -- which is additional and separate from the definition $\tilde{\mathbf{u}}$ in (A9). This alteration of $\tilde{\mathbf{u}}$ on the boundary implicitly contributes the appropriate terms that account for inhomogeneous boundary-normal pressure gradients $\hat{n} \cdot \nabla p_{nh}^{m+1} \ne 0$ to the right-hand-side of (A12) during the computation of $\nabla \cdot \tilde{\mathbf{u}}$.

A preconditioned conjugate gradient iteration may be used on non-rectilinear grids, including complex domains. For domains that immerse an irregular boundary within a RectilinearGrid, we have implemented an efficient, rapidly-converging preconditioner that leverages the FFT-based solver with masking applied to immersed cells. The FFT-based preconditioner for solving the Poisson equation in irregular domains will be described in a forthcoming paper.

### A3 Time discretization of the HydrostaticFreeSurfaceModel

The HydrostaticFreeSurfaceModel uses a linear free surface formulation paired with a geopotential vertical coordinate that may be integrated in time using either a fully ExplicitFreeSurface, an ImplicitFreeSurface utilizing a two-dimensional elliptical solve, or a SplitExplicitFreeSurface. The latter free surface solver can also be used to solve the primitive equations with a non-linear free surface formulation and a free-surface following vertical coordinate (the z-star vertical coordinate, Adcroft & Campin, 2004). For brevity, we describe here only the SplitExplicitFreeSurface, which is the most generally useful method. The SplitExplicitFreeSurface substeps the depth-integrated or "barotropic" horizontal velocity $\mathbf{U}_h$ along with the free surface displacement $\eta$ using a short time step while and the depth-dependent, "baroclinic" velocities, along with tracers, are relatively stationary.

The barotropic horizontal transport $\mathbf{U}_h$ is defined

$$\mathbf{U}_h \stackrel{\text{def}}{=} \int_{-H}^{\eta} \mathbf{u}_h \, dz \,, \tag{A13}$$

where $\mathbf{u}_h = (u, v)$ is the total horizontal velocity and $H$ is the depth of the fluid. Similarly integrating the horizontal momentum equations (17) from $z = -H$ to $z = \eta$ yields an evolution equation for $\mathbf{U}_h$,

$$\partial_t \mathbf{U}_h = -g(H + \eta)\nabla_h \eta + \int_{-H}^{\eta} \mathbf{G}_{u_h} \, dz \,, \tag{A14}$$

where $\mathbf{G}_{u_h}$ includes all the tendency terms that evolve "slowly" compared to the barotropic mode:

$$\mathbf{G}_{u_h} = -(\mathbf{u} \cdot \nabla)\mathbf{u}_h - \mathbf{f} \times \mathbf{u} - \nabla \cdot \boldsymbol{\tau} + \mathbf{F}_h \,. \tag{A15}$$

The evolution equation for the free surface is obtained by integrating the continuity equation (15) in $z$ to obtain $\nabla \cdot \mathbf{U}_h = -w|_{z=\eta}$, and inserting this into (16) to find

$$\partial_t \eta = -\nabla_h \cdot \mathbf{U}_h \,. \tag{A16}$$

The pair of equations (A14) and (A16) characterize the evolution of the barotropic mode, which involves faster time-scales than the baroclinic mode evolution described by equations (17). To resolve both modes, we use a split-explicit algorithm where the barotropic mode is advanced in time using a smaller time-step than the one used for three-dimensional baroclinic variables. In particular, a predictor three-dimensional velocity is evolved without accounting for the barotropic mode evolution, using the QAB2 scheme described by section A3. We denote this "predictor" velocity, again, with a tilde as done in section A2.

$$(1 - \Delta t\, \partial_z \nu_z^m \partial_z)\, \tilde{\mathbf{u}}_h - \mathbf{u}_h^m \approx \int_{t^m}^{t^{m+1}} \mathbf{G}_{u_h} \, dt \,. \tag{A17}$$

We then compute the barotropic mode evolution by sub-stepping $M$ times the barotropic equations using a forward-backward time-stepping scheme and a time-step $\Delta\tau = \Delta t / N$,

$$\eta^{n+1} - \eta^n = -\Delta\tau\, \nabla_h \cdot \mathbf{U}_h^n \,, \tag{A18}$$

$$\mathbf{U}_h^{n+1} - \mathbf{U}_h^n = -\Delta\tau \left[ g(H + \eta)\nabla_h \eta^{n+1} - \frac{1}{\Delta t} \int_{-H}^{\eta} \int_{t^m}^{t^{m+1}} \mathbf{G}_{u_h} \, dt \, dz \right] \,. \tag{A19}$$

The slow tendency terms are frozen in time during substepping. The barotropic quantities are averaged within the sub-stepping with

$$\bar{\mathbf{U}}_h = \sum_{n=1}^{M} a_n \mathbf{U}_h^n \quad \text{and} \quad \bar{\eta} = \sum_{n=1}^{M} a_n \eta^n \,, \tag{A20}$$

where $M$ is the number of substeps per baroclinic step and $a_n$ are the weights calculated from the provided averaging kernel. The default choice of averaging kernel is the minimal dispersion filters developed by Shchepetkin and McWilliams (2005). The number of substeps $M$ is calculated to center the averaging kernel at $t^{m+1}$. As a result, the barotropic subcycling overshoots the baroclinic step, i.e. $M > N$ with a maximum of $M = 2N$. Finally, the barotropic mode is reconciled to the baroclinic mode with a correction step

$$\mathbf{u}_h^{m+1} = \tilde{\mathbf{u}}_h + \frac{1}{H + \eta}\left(\bar{\mathbf{U}}_h - \int_{-H}^{\eta} \tilde{\mathbf{u}}_h \, dz\right) \,. \tag{A21}$$

The barotropic variables are then reinitialized for evolution in the next barotropic mode evolution using the time-averaged $\bar{\eta}$ and $\bar{\mathbf{U}}_h$.

## Appendix B: Finite volume spatial discretization

Oceananigans uses a finite volume method in which fields are represented discretely by their average value over small local regions or "finite volumes" of the domain. Listing 13 discretizes $c = e^{xy}$ on three different grids that cover the unit square. At the finest resolution, each cell-averaged value $c_{ij}^{\text{fine}}$ is computed approximately using `set!` to evaluate $e^{xy}$ at the center of each finite volume, where $i, j$ denote the $x$ and $y$ indices of the finite volumes. At medium and coarse resolution, the $c_{ij}^{\text{medium}}$ and $c_{ij}^{\text{coarse}}$ are computed by averaging or "regridding" fields discretized at a higher resolution. This computation produces three fields with identical integrals over the unit square. For example, integrals are computed exactly by summing discrete fields over all cells,

<!-- [Equation B1] -->
$$\int c \, dx \, dy = \sum_{i,j}^{1024,1024} V_{ij}^{\text{fine}} c_{ij}^{\text{fine}} = \sum_{i,j}^{16,16} V_{ij}^{\text{medium}} c_{ij}^{\text{medium}} = \sum_{i,j}^{4,4} V_{ij}^{\text{coarse}} c_{ij}^{\text{coarse}} \,, \tag{B1}$$

where $V_{ij}$ is the "volume" of the cell with indices $i, j$ (more accurately an "area" in this two-dimensional situation). Figure B1 visualizes the three fields.

**Listing 13:** Finite volume discretization of $e^{xy}$ on three grids over the unit square. The fields are visualized in figure B1. The meaning of the "Center" in "CenterField" is discussed below.

```julia
topology = (Bounded, Bounded, Flat)
x = y = (0, 1)
c(x, y) = exp(x) * y

fine_grid = RectilinearGrid(size=(1024, 1024); x, y, topology)
c_fine = CenterField(fine_grid)
set!(c_fine, c)

medium_grid = RectilinearGrid(size=(16, 16); x, y, topology)
c_medium = CenterField(medium_grid)
regrid!(c_medium, c_fine)

coarse_grid = RectilinearGrid(size=(4, 4); x, y, topology)
c_coarse = CenterField(coarse_grid)
regrid!(c_coarse, c_medium)
```

[Figure B1: Finite volume discretization of e^xy on the unit square at three different resolutions.]

The discrete calculus and arithmetic operations required to solve the governing equations of the NonhydrostaticModel and HydrostaticFreeSurfaceModel use the system of "staggered grids" described by Arakawa and Lamb (1977). Both models use "C-grid" staggering, where cells for tracers, pressure, and the divergence of the velocity field $\nabla \cdot \mathbf{u}$ are co-located, and cells for velocity components $\mathbf{u} = (u, v, w)$ are staggered by half a cell width in the x-, y-, and z-direction, respectively. Listing 14 illustrates grid construction and notation for a one-dimensional staggered grid with unevenly-spaced cells. Figure B2 visualizes 2- and 3-dimensional staggered grids, indicating the location of certain variables.

**Listing 14:** A one-dimensional staggered grid.

```julia
using Oceananigans

grid = RectilinearGrid(topology=(Bounded, Flat, Flat), size=4, x=[0, 0.2, 0.3, 0.7, 1])

u = Field{Face, Center, Center}(grid)
c = Field{Center, Center, Center}(grid)

xnodes(u)    # [0.0, 0.2, 0.3, 0.7, 1.0]
xnodes(c)    # [0.1, 0.25, 0.5, 0.85]
location(dx(c)) # (Face, Center, Center)
```

[Figure B2: Locations of cell centers and interfaces on a two-dimensional (a) and three-dimensional (b) staggered grid. In (a), the red and blue shaded regions highlight the volumes in the dual u-grid and v-grid, located at (Face, Center, Center) and (Center, Face, Center), respectively. In (b), the shaded regions highlight the facial areas used in the fluxes computations, denoted with Ax, Ay, and Az.]

### B1 A system of composable operators

A convention for indexing is associated with staggered locations. Face indices are "left" of cell indices. This means that difference operators acting on fields at cells differ from those that act on face fields. To illustrate this we introduce Oceananigans-like difference operators,

```julia
delta_x_fcc(i, j, k, grid, c) = c[i, j, k] - c[i-1, j, k]
delta_x_ccc(i, j, k, grid, u) = u[i+1, j, k] - u[i, j, k]
```

where superscripts denote the location of the result of the operation. For example, the difference $\delta_x^{fcc}$ acts on fields located at $ccc$ (meaning cell Center in the x, y, and z directions respectively). Complementary to the difference operators are reconstruction or "interpolation" operators,

```julia
interp_x_fcc(i, j, k, grid, c) = (c[i, j, k] + c[i-1, j, k]) / 2
interp_x_ccc(i, j, k, grid, u) = (u[i+1, j, k] + u[i, j, k]) / 2
```

The prefix arguments `i, j, k, grid` are more than convention: the prefix enables a system for composing operators. For example, defining

```julia
delta_x_fcc(i, j, k, grid, f::Function, args...) =
    f(i, j, k, grid, args...) - f(i-1, j, k, grid, args...)

delta_x_ccc(i, j, k, grid, f::Function, args...) =
    f(i+1, j, k, grid, args...) - f(i, j, k, grid, args...)
```

leads to a concise definition of the second-difference operator:

```julia
delta2_x_ccc(i, j, k, grid, f::Function, a...) = delta_x_ccc(i, j, k, grid, delta_x_fcc, f, a...)
```

Operator composition is used throughout Oceananigans source code to implement stencil operations.

### B2 Tracer flux divergences, advection schemes, and reconstruction

The divergence of a tracer flux $\mathbf{J} = J_x\,\hat{x} + J_y\,\hat{y} + J_z\,\hat{z}$ is discretized conservatively by the finite volume method via

$$\nabla \cdot \mathbf{J} \approx \frac{1}{V_c} \left[ \delta_x\left(\underbrace{A_x J_x}_{fcc}\right) + \delta_y\left(\underbrace{A_y J_y}_{cfc}\right) + \delta_z\left(\underbrace{A_z J_z}_{ccf}\right) \right] \,, \tag{B2}$$

where $\delta_x, \delta_y, \delta_z$ are difference operators in $x, y, z$, $V_c$ denotes the volume of the tracer cells, $A_x$, $A_y$, and $A_z$ denote the areas of the tracer cell faces with surface normals $\hat{x}$, $\hat{y}$, and $\hat{z}$, respectively. Equation (B2) indicates the location of each flux component: fluxes into tracers cell at $ccc$ are computed at the cell faces located at $fcc$, $cfc$, and $ccf$.

The advective tracer flux in (9) is written in "conservative form" using incompressibility (2), and then discretized similarly to (B2) to form

$$\mathbf{u} \cdot \nabla c = \nabla \cdot (\mathbf{u}\,c) \approx \frac{1}{V_c} \left[ \delta_x\left(A_x u\,\langle c \rangle_x\right) + \delta_y\left(A_y v\,\langle c \rangle_y\right) + \delta_z\left(A_z w\,\langle c \rangle_z\right) \right] \,, \tag{B3}$$

where $\langle c \rangle_x$ denotes a reconstruction of $c$ in the x-direction from its native location $ccc$ to the tracer cell interface at $fcc$; $\langle c \rangle_y$ and $\langle c \rangle_z$ in (B3) are defined similarly.

The advective fluxes $\mathbf{u}\,c$ must be computed on interfaces between tracer cells, where the approximate value of $c$ must be reconstructed. (Velocity components like $u$ must also be reconstructed on interfaces. Within the C-grid framework, we approximate $u$ on tracer cell interfaces directly using the values $u_{ijk}$, which represent $u$ averaged over a region encompassing the interface.) The simplest kind of reconstruction is `Centered(order=2)`, which is equivalent to taking the average between adjacent cells,

$$\langle c \rangle_i = \frac{1}{2}(c_i + c_{i-1}) \,, \tag{B4}$$

where $\langle c \rangle_i$ denotes the centered reconstruction of $c$ on the interface at $x = x_{i-1/2}$. Also in (B4) the $j, k$ indices are implied and we have suppressed the direction $x$ to lighten the notation. Reconstructions stencils for `Centered(order=N)` are automatically generated for even $N$ up to $N_{\max} = 12$, where $N_{\max}$ is an adjustable parameter in the source code. All subsequent reconstructions are described in the x-direction only.

Centered schemes should be used when explicit dissipation justified by a physical rationale dominates at the grid scale. In scenarios where dissipation is needed solely for artificial reasons, we find applications for UpwindBiased schemes, which use an odd-order stencil biased against the direction of flow. For example, `UpwindBiased(order=1)` and `UpwindBiased(order=3)` schemes are written

<!-- [Equation B5] -->
$$u[c]_x^1 = \begin{cases} u\,c_{i-1} & \text{if } u > 0 \\ u\,c_i & \text{if } u < 0 \end{cases} \quad \text{and} \quad u[c]_x^3 = \begin{cases} u\,\frac{1}{6}(-c_{i-2} + 5c_{i-1} + 2c_i) & \text{if } u > 0 \\ u\,\frac{1}{6}(2c_{i-1} + 5c_i - c_{i+1}) & \text{if } u < 0 \end{cases} \tag{B5}$$

where $[c]_x^N$ denotes $N$th-order upwind reconstruction in the x-direction. (Note that $u[c]_x^N = 0$ if $u = 0$.)

The compact form of equations (B5) demonstrates how upwind schemes introduce variance dissipation through numerical discretization. In particular, an `UpwindBiased(order=1)` reconstruction can be rewritten as a sum of a `Centered(order=2)` discrete advective flux and a discrete diffusive flux

$$u[c]_x^1 = u\frac{c_i + c_{i-1}}{2} - \kappa_1 \frac{c_i - c_{i-1}}{\Delta x} \,, \quad \text{where} \quad \kappa_1 = \frac{|u|\Delta x}{2} \,. \tag{B6}$$

Reordering the `UpwindBiased(order=3)` advective flux in the same manner recovers a sum of a `Centered(order=4)` advective flux and a 4th-order hyperdiffusive flux, equivalent to a finite volume approximation of

$$uc + \kappa_3 \frac{\partial^3 c}{\partial x^3} \,, \quad \text{where} \quad \kappa_3 = \frac{|u|\Delta x^3}{12} \,. \tag{B7}$$

UpwindBiased reconstruction can be always reordered to expose a sum of Centered reconstruction and a high-order diffusive flux with a velocity-dependent diffusivity. The diffusive operator associated with `UpwindBiased(order=1)` and `UpwindBiased(order=3)` is enough to offset the dispersive errors of the Centered component and, therefore, eliminate the artificial explicit diffusion needed for stability purposes. However, this approach does not scale to high order since the diffusive operator associated with a high order UpwindBiased scheme (5th, 7th, and so on), becomes quickly insufficient to eliminate spurious errors associated with the Centered component (Godunov, 1959).

The inability to achieve high order and, therefore, low dissipation motivated the implementation of Weighted, Essentially Non-Oscillatory (WENO) reconstruction (C. Shu, 1997; C.-W. Shu, 2009). WENO is a non-linear reconstruction scheme that combines a set of odd-order linear reconstructions obtained by stencils that are shifted by a value $s$ relative to the canonical UpwindBiased stencil, using a weighting scheme for each stencil that depends on the smoothness of the reconstructed field $c$. Since the constituent stencils are lower-order than the WENO order, this strategy yields a scheme whose order of accuracy adapts depending on the smoothness of the reconstructed field. In smooth regions high-order is retained, while the order quickly decreases in the presence of noisy regions, decreasing the order of the associated diffusive operator. WENO proves especially useful for high-resolution, turbulence-resolving simulations (either at meter or planetary scales) without requiring any additional explicit artificial dissipation (Pressel et al., 2017; Silvestri et al., 2024).

To illustrate how WENO works we consider a fifth-order WENO scheme for $u > 0$,

$$\{c\}^5 = \gamma_0 [c]^{3,0} + \gamma_1 [c]^{3,1} + \gamma_2 [c]^{3,2} \,, \tag{B8}$$

where the notation $[c]^{3,s}$ denotes an UpwindBiased stencil shifted by $s$ indices, such that $[c]^3 \stackrel{\text{def}}{=} [c]^{3,0}$. The shifted upwind stencils $[c]_i^{N,s}$ evaluated at index $i$ are defined

$$[c]_i^{3,s} = \frac{1}{6} \begin{cases} -c_{i-1} + 5c_i + 2c_{i+1} & \text{for } s = -1 \\ 2c_{i-2} + 5c_{i-1} - c_i & \text{for } s = 0 \\ 2c_{i-3} - 7c_{i-2} + 11c_{i-1} & \text{for } s = 1 \end{cases} \tag{B9}$$

The weights $\gamma_s(c)$ are determined by a smoothness metric that produces $\{c\}^5 \approx [c]^5$ when $c$ is smooth, but limits to the more diffusive $\{c\}^5 \approx [c]^3$ when $c$ changes abruptly. Thus WENO adaptively introduces dissipation as needed based on the smoothness of $c$, yielding stable simulations with a high effective resolution that require no artificial dissipation. WENO can alternatively be interpreted as adding an implicit hyperviscosity that adapts from low- to high-order depending on the local nature of the solution. To compute the weights $\gamma_s(c)$, we use the WENO-Z formulation (Balsara & Shu, 2000).

The properties of Centered, UpwindBiased, and WENO reconstruction are investigated by listing 15, which simulates the advection of a top hat tracer distribution. The results are plotted in figure B3.

**Listing 15:** A script that advects a top hat tracer profile in one-dimension with a constant prescribed velocity. We use halo=6 to accommodate schemes up to WENO(order=11).

```julia
using Oceananigans

grid = RectilinearGrid(size=128; x=(-4, 8), halo=6, topology=(Periodic, Flat, Flat))
advection = WENO(order=9) # Centered(order=2), UpwindBiased(order=3)
velocities = PrescribedVelocityFields(u=1)
model = HydrostaticFreeSurfaceModel(; grid, velocities, advection, tracers=:c)

top_hat(x) = abs(x) > 1 ? 0 : 1
set!(model, c = top_hat)

simulation = Simulation(model, dt=1/grid.Nx, stop_time=4)
run!(simulation)
```

[Figure B3: Advection of a top hat tracer distribution in one-dimension using various advection schemes.]

#### B2.1 Discretization of momentum advection

The discretization of momentum advection with a flux form similar to (B3) is more complex than the tracer case because both the advecting velocity and advected velocity require reconstruction. We use the method described by Ghosh and Baeder (2012) and Pressel et al. (2015), wherein advecting velocities are constructed with a high-order Centered scheme when the advected velocity component is reconstructed with a high-order UpwindBiased or WENO scheme. We have also developed a novel WENO-based method for discretizing momentum advection in the rotational or "vector invariant" form especially appropriate for representing mesoscale and submesoscale turbulent advection on curvilinear grids (Silvestri et al., 2024).

## Appendix C: Parallelization

Oceananigans supports distributed computations with slab and pencil domain decomposition. The interior domain is extended using "halo" or "ghost" cells that hold the results of interprocessor boundaries. "halo" cells are updated before the computation of tendencies through asynchronous send / receive operations using the message passing interface (MPI) Julia library (Byrne et al., 2021). For a detailed description of the parallelization strategy of the HydrostaticFreeSurfaceModel; see Silvestri et al. (2025). The NonhydrostaticModel implements the same overlap of communication and computation for halo exchange before the calculation of tendencies. For the FFT-based three-dimensional pressure solver, we implement a transpose algorithm that switches between x-local, y-local, and z-local configurations to compute efficiently the discrete transforms. The transpose algorithm for the distributed FFT solver is shown in figure C1.

[Figure C1: A schematic showing the distributed Poisson solver procedure with a pencil parallelization that divides the domain in two ranks in both x and y. The schematic highlights the data layout in the ranks during each operation.]

## Appendix D: Table of numerical examples

Table D1 provides a list of the numerical experiments in this paper.

| Description | Code | Visualization |
|---|---|---|
| 2D turbulence using WENO(order=9) advection | listing 1 | fig. 1 |
| 2D turbulence with moving tracer source | listing 2 | fig. 1 |
| DNS and LES of flow around a cylinder at various Re | listing 4 | fig. 2 |
| DNS of cabbeling in freshwater | listing 5 | fig. 3 |
| LES of the Eady problem with WENO(order=9) | listing 6 | fig. 4 |
| Tidally-oscillating flow past Three Tree Point | listing 7 | fig. 5 |
| Internal waves generated by tidal forcing over bathymetry | listing 8 | fig. 6 |
| Comparison of vertical mixing parameterizations | listing 9 | fig. 7 |
| Near-global baroclinic instability on three grids | listing 10 | fig. 8 |
| Realistic ocean simulation with ClimaOcean | listing 11 | fig. 9 |
| Tracer diffusion with various time discretizations | listing 12 | fig. A1 |
| Visualization of the finite volume discretization | listing 13 | fig. B1 |
| One-dimensional advection of a top-hat tracer profile | listing 15 | fig. B3 |

## Open Research Section

Oceananigans is available at the GitHub repository github.com/CliMA/Oceananigans.jl and ClimaOcean is available at github.com/CliMA/ClimaOcean.jl. Oceananigans documentation lives at clima.github.io/OceananigansDocumentation. All the scripts that reproduce the simulations and figures in this paper are available at the GitHub repository github.com/glwagner/OceananigansPaper. Visualizations were made using Makie.jl (Danisch & Krumbiegel, 2021).

## Acknowledgments

This project is supported by Schmidt Sciences, LLC and by the National Science Foundation grant AGS-1835576. N.C.C. is additionally supported by the Australian Research Council under the Center of Excellence for the Weather of the 21st Century CE230100012 and the Discovery Project DP240101274.

## References

Abkar, M., Bae, H. J., & Moin, P. (2016). Minimum-dissipation scalar transport model for large-eddy simulation of turbulent flows. *Physical Review Fluids*, 1(4), 041701.

Adcroft, A., Anderson, W., Balaji, V., Blanton, C., Bushuk, M., Dufour, C. O., ... Zhang, R. (2019). The GFDL global ocean and sea ice model OM4.0: Model description and simulation features. *Journal of Advances in Modeling Earth Systems*, 11(10), 3167-3211. doi: 10.1029/2019MS001726

Adcroft, A., & Campin, J.-M. (2004). Rescaled height coordinates for accurate representation of free-surface flows in ocean circulation models. *Ocean Modelling*, 7(3-4), 269-284.

Adcroft, A., Campin, J.-M., Hill, C., & Marshall, J. (2004). Implementation of an atmosphere-ocean general circulation model on the expanded spherical cube. *Monthly Weather Review*, 132(12), 2845-2863.

Adcroft, A., Hill, C., & Marshall, J. (1997). Representation of topography by shaved cells in a height coordinate ocean model. *Monthly Weather Review*, 125(9), 2293-2315.

Alnaes, M., Blechta, J., Hake, J., Johansson, A., Kehlet, B., Logg, A., ... Wells, G. N. (2015). The FEniCS project version 1.5. *Archive of Numerical Software*, 3(100).

Arakawa, A., & Lamb, V. R. (1977). Computational design of the basic dynamical processes of the UCLA general circulation model.

Auclair, F., Benshila, R., Bordois, L., Boutet, M., Bremond, M., Caillaud, M., ... Zribi, A. (2025). Coastal and Regional Ocean COmmunity model. *Zenodo*. doi: 10.5281/zenodo.7415055

Balsara, D., & Shu, C. (2000). Monotonicity preserving weighted essentially non-oscillatory schemes with increasingly high order of accuracy. *Journal of Computational Physics*, 160(2), 405-452. doi: 10.1006/jcph.2000.6443

Barnes, A. J., Constantinou, N. C., Gibson, A. H., Kiss, A. E., Chapman, C., Reily, J., ... Yang, L. (2024). regional-mom6: A Python package for automatic generation of regional configurations for the Modular Ocean Model 6. *Journal of Open Source Software*, 9(100), 6857. doi: 10.21105/joss.06857

Besard, T., Foket, C., & De Sutter, B. (2018). Effective extensible programming: unleashing Julia on GPUs. *IEEE Transactions on Parallel and Distributed Systems*, 30(4), 827-841.

Bezanson, J., Edelman, A., Karpinski, S., & Shah, V. B. (2017). Julia: A fresh approach to numerical computing. *SIAM review*, 59(1), 65-98. doi: 10.1137/141000671

Bisits, J. I., Zika, J. D., & Sohail, T. (2025). Cabbeling as a catalyst and driver of turbulent mixing. *Journal of Fluid Mechanics*, 1011, A17. doi: 10.1017/jfm.2025.349

Boccaletti, G., Ferrari, R., & Fox-Kemper, B. (2007). Mixed layer instabilities and restratification. *Journal of Physical Oceanography*, 37(9), 2228-2250. doi: 10.1175/JPO3101.1

Bou-Zeid, E., Meneveau, C., & Parlange, M. (2005). A scale-dependent Lagrangian dynamic model for large eddy simulation of complex turbulent flows. *Physics of fluids*, 17(2).

Bradbury, J., Frostig, R., Hawkins, P., Johnson, M. J., Leary, C., Maclaurin, D., ... Zhang, Q. (2018). JAX: composable transformations of Python+NumPy programs. Retrieved from http://github.com/jax-ml/jax

Bryan, K. (1969). A numerical method for the study of the circulation of the world ocean. *Journal of Computational Physics*, 135(2), 154-169.

Burns, K. J., Vasil, G. M., Oishi, J. S., Lecoanet, D., & Brown, B. P. (2020). Dedalus: A flexible framework for numerical simulations with spectral methods. *Physical Review Research*, 2(2), 023068. doi: 10.1103/PhysRevResearch.2.023068

Byrne, S., Wilcox, L. C., & Churavy, V. (2021). MPI.jl: Julia bindings for the Message Passing Interface. In *Proceedings of the JuliaCon Conferences* (Vol. 1, p. 68). doi: 10.21105/jcon.00068

Callies, J., & Ferrari, R. (2018). Baroclinic instability in the presence of convection. *Journal of Physical Oceanography*, 48(1), 45-60. doi: 10.1175/JPO-D-17-0028.1

Chassignet, E. P., & Xu, X. (2017). Impact of horizontal resolution (1/12 to 1/50) on Gulf Stream separation, penetration, and variability. *Journal of Physical Oceanography*, 47(8), 1999-2021.

Chassignet, E. P., & Xu, X. (2021). On the importance of high-resolution in large-scale ocean models. *Advances in Atmospheric Sciences*, 38, 1621-1634.

Chor, T., Constantinou, N. C., Bisits, J. I., Wagner, G. L., Ramadhan, A. R., Zheng, Z., & Whitley, V. (2025). Oceanostics.jl. *Zenodo*. doi: 10.5281/zenodo.8280754

Churavy, V. (2024). KernelAbstractions.jl. *Zenodo*. doi: 10.5281/zenodo.13773520

Cox, M. D. (1984). A primitive equation, 3-dimensional model of the ocean (Tech. Rep. No. 1). Princeton, NJ: NOAA Geophysical Fluid Dynamics Laboratory.

Craik, A. D., & Leibovich, S. (1976). A rational model for Langmuir circulations. *Journal of Fluid Mechanics*, 73(3), 401-426.

Danilov, S., Sidorenko, D., Wang, Q., & Jung, T. (2017). The finite-volume sea ice-ocean model (FESOM2). *Geoscientific Model Development*, 10(2), 765-789.

Danisch, S., & Krumbiegel, J. (2021). Makie.jl: Flexible high-performance data visualization for Julia. *Journal of Open Source Software*, 6(65), 3349. doi: 10.21105/joss.03349

Dansereau, V., Weiss, J., Saramito, P., & Lattes, P. (2016). A Maxwell elasto-brittle rheology for sea ice modelling. *The Cryosphere*, 10(3), 1339-1359.

Dheeshjith, S., Subel, A., Adcroft, A., Busecke, J., Fernandez-Granda, C., Gupta, S., & Zanna, L. (2025). Samudra: An AI global ocean emulator for climate. *Geophysical Research Letters*, 52(10), e2024GL114318.

Dong, J., Fox-Kemper, B., Wenegrat, J. O., Bodner, A. S., Yu, X., Belcher, S., & Dong, C. (2024). Submesoscales are a significant turbulence source in global ocean surface boundary layer. *Nature Communications*, 15(1), 9566.

Eady, E. T. (1949). Long waves and cyclone waves. *Tellus*, 1(3), 33-52.

Edson, J. B., Jampana, V., Weller, R. A., Bigorre, S. P., Plueddemann, A. J., Fairall, C. W., ... Hersbach, H. (2014). On the exchange of momentum over the open ocean. *Journal of Physical Oceanography*, 44(9), 1589.

Forget, G., Campin, J.-M., Heimbach, P., Hill, C. N., Ponte, R. M., & Wunsch, C. (2015). ECCO version 4: An integrated framework for non-linear inverse modeling and global ocean state estimation. *Geoscientific Model Development*, 8(10), 3071-3104.

Garrett, C., & Kunze, E. (2007). Internal tide generation in the deep ocean. *Annu. Rev. Fluid Mech.*, 39(1), 57-87.

Gent, P. R., & Mcwilliams, J. C. (1990). Isopycnal mixing in ocean circulation models. *Journal of Physical Oceanography*, 20(1), 150-155. doi: 10.1175/1520-0485(1990)020<0150:IMIOCM>2.0.CO;2

Ghosh, D., & Baeder, J. D. (2012). High-order accurate incompressible Navier-Stokes algorithm for vortex-ring interactions with solid wall. *AIAA journal*, 50(11), 2408-2422.

Godunov, S. K. (1959). A difference scheme for numerical solution of discontinuous solution of hydrodynamic equations. *Matematicheskii Sbornik*, 47, 271-306.

Good, S. A., Martin, M. J., & Rayner, N. A. (2013). EN4: Quality controlled ocean temperature and salinity profiles and monthly objective analyses with uncertainty estimates. *Journal of Geophysical Research: Oceans*, 118(12), 6704-6716.

Grachev, A. A., Andreas, E. L., Fairall, C. W., Guest, P. S., & Persson, P. O. G. (2007). SHEBA flux-profile relationships in the stable atmospheric boundary layer. *Boundary-layer meteorology*, 124, 315-333.

Griffies, S. M. (1998). The Gent-McWilliams skew flux. *Journal of Physical Oceanography*, 28(5), 831-841. doi: 10.1175/1520-0485(1998)028<0831:TGMSF>2.0.CO;2

Griffies, S. M., Adcroft, A., & Hallberg, R. W. (2020). A primer on the vertical Lagrangian-remap method in ocean models based on finite volume generalized vertical coordinates. *Journal of Advances in Modeling Earth Systems*, 12(10), e2019MS001954.

Griffies, S. M., Pacanowski, R. C., & Hallberg, R. W. (2000). Spurious diapycnal mixing associated with advection in a z-coordinate ocean model. *Monthly Weather Review*, 128(3), 538-564.

Griffies, S. M., Stouffer, R. J., Adcroft, A. J., Bryan, K., Dixon, K. W., Hallberg, R., ... Rosati, A. (2015). A historical introduction to MOM.

Hafner, D., Nuterman, R., & Jochum, M. (2021). Fast, cheap, and turbulent -- global ocean modeling with GPU acceleration in Python. *Journal of Advances in Modeling Earth Systems*, 13(12), e2021MS002717.

Halliwell, G. R. (2004). Evaluation of vertical coordinate and vertical mixing algorithms in the HYbrid-Coordinate Ocean Model (HYCOM). *Ocean Modelling*, 7(3-4), 285-322.

Harcourt, R. R. (2015). An improved second-moment closure model of Langmuir turbulence. *Journal of Physical Oceanography*, 45(1), 84-103.

Harlow, F. H., Welch, J. E., et al. (1965). Numerical calculation of time-dependent viscous incompressible flow of fluid with free surface. *Physics of fluids*, 8(12), 2182.

Held, I. M. (2005). The gap between simulation and understanding in climate modeling. *Bulletin of the American Meteorological Society*, 86(11), 1609-1614.

Hell, M., Chapron, B., & Fox-Kemper, B. (2024). A particle-in-cell wave model for efficient sea-state and swell estimates in Earth system models -- PiCLES. *Authorea Preprints*.

Huang, N. E. (1979). On surface drift currents in the ocean. *Journal of Fluid Mechanics*, 91(1), 191-208.

Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9(3), 90-95. doi: 10.1109/MCSE.2007.55

Jansen, M. F., Adcroft, A., Khani, S., & Kong, H. (2019). Toward an energetically consistent, resolution aware parameterization of ocean mesoscale eddies. *Journal of Advances in Modeling Earth Systems*, 11(8), 2844-2860.

Jullien, S., Le Corre, M., Raynaud, S., Gula, J., Sepulveda, A., Le Gentil, S., ... Collin, J. (2025). CROCO Pytools. *Zenodo*. doi: 10.5281/zenodo.15064146

Karna, T., Kramer, S. C., Mitchell, L., Ham, D. A., Piggott, M. D., & Baptista, A. M. (2018). Thetis coastal ocean model: discontinuous Galerkin discretization for the three-dimensional hydrostatic equations. *Geoscientific Model Development*, 11(11), 4359-4382.

Khairoutdinov, M. F., Blossey, P. N., & Bretherton, C. S. (2022). Global system for atmospheric modeling: Model description and preliminary results. *Journal of Advances in Modeling Earth Systems*, 14(6), e2021MS002968.

Kimmritz, M., Danilov, S., & Losch, M. (2016). The adaptive EVP method for solving the sea ice momentum equation. *Ocean Modelling*, 101, 59-67.

Kiss, A. E., Hogg, A. M., Hannah, N., Boeira Dias, F., Brassington, G. B., Chamberlain, M. A., ... others (2020). ACCESS-OM2 v1.0: a global ocean-sea ice model at three resolutions. *Geoscientific Model Development*, 13(2), 401-442.

Klower, M., Gelbrecht, M., Hotta, D., Willmert, J., Silvestri, S., Wagner, G. L., ... others (2024). SpeedyWeather.jl: Reinventing atmospheric general circulation models towards interactivity and extensibility. *Journal of Open Source Software*, 9(98), 6323. doi: 10.21105/joss.06323

Kochkov, D., Smith, J. A., Alieva, A., Wang, Q., Brenner, M. P., & Hoyer, S. (2021). Machine learning-accelerated computational fluid dynamics. *Proceedings of the National Academy of Sciences*, 118(21), e2101784118.

Kochkov, D., Yuval, J., Langmore, I., Norgaard, P., Smith, J., Mooers, G., ... Hoyer, S. (2024). Neural general circulation models for weather and climate. *Nature*, 632, 1060-1066. doi: 10.1038/s41586-024-07744-y

Korn, P., Bruggemann, N., Jungclaus, J. H., Lorenz, S., Gutjahr, O., Haak, H., ... others (2022). ICON-O: The ocean component of the ICON Earth system model -- Global simulation characteristics and local telescoping capability. *Journal of Advances in Modeling Earth Systems*, 14(10), e2021MS002952.

Lam, R., Sanchez-Gonzalez, A., Willson, M., Wirnsberger, P., Fortunato, M., Alet, F., ... others (2023). Learning skillful medium-range global weather forecasting. *Science*, 382(6677), 1416-1421.

Leclair, M., & Madec, G. (2011). z-Coordinate, an arbitrary Lagrangian-Eulerian coordinate separating high and low frequency motions. *Ocean Modelling*, 37(3-4), 139-152.

Lecoanet, D., McCourt, M., Quataert, E., Burns, K. J., Vasil, G. M., Oishi, J. S., ... O'Leary, R. M. (2016). A validated non-linear Kelvin-Helmholtz benchmark for numerical hydrodynamics. *Monthly Notices of the Royal Astronomical Society*, 455(4), 4274-4288.

Legay, A., Deremble, B., & Burchard, H. (2025). Derivation and implementation of a non-local term to improve the oceanic convection representation within the k-epsilon parameterization. *Journal of Advances in Modeling Earth Systems*, 17(1), e2024MS004243. doi: 10.1029/2024MS004243

Lilly, D. K. (1983). Stratified turbulence and the mesoscale variability of the atmosphere. *Journal of the Atmospheric Sciences*, 40(3), 749-761. doi: 10.1175/1520-0469(1983)040<0749:STATMV>2.0.CO;2

Mak, J., Maddison, J. R., Marshall, D. P., & Munday, D. R. (2018). Implementation of a geometrically informed and energetically constrained mesoscale eddy parameterization in an ocean circulation model. *Journal of Physical Oceanography*, 48(10), 2363-2382.

Marshall, J., Adcroft, A., Hill, C., Perelman, L., & Heisey, C. (1997). A finite-volume, incompressible Navier Stokes model for studies of the ocean on parallel computers. *Journal of Geophysical Research: Oceans*, 102(C3), 5753-5766.

McDougall, T. J., & Barker, P. M. (2011). Getting started with TEOS-10 and the Gibbs Seawater (GSW) oceanographic toolbox. *SCOR/IAPSO WG*, 127(532), 1-28.

Molemaker, M. J., McWilliams, J. C., & Capet, X. (2010). Balanced and unbalanced routes to dissipation in an equilibrated Eady flow. *Journal of Fluid Mechanics*, 654, 35-63.

Monin, A. A., & Obukhov, A. M. (1954). Basic laws of turbulent mixing in the surface layer of the atmosphere. *Contrib. Geophys. Inst. Acad. Sci. USSR*, 151(163), e187.

Moses, W. S., Churavy, V., Paehler, L., Huckelheim, J., Narayanan, S. H. K., Schanen, M., & Doerfert, J. (2021). Reverse-mode automatic differentiation and optimization of GPU kernels via Enzyme. In *Proceedings of the international conference for high performance computing, networking, storage and analysis* (pp. 1-16).

Mrozowska, M. A., Avery, J., Stoustrup, A., Nuterman, R., Johnsen, C.-J., Thormann, A., & Jochum, M. (2025). Bayesian optimization with GPU acceleration for ocean models. *Journal of Geophysical Research: Machine Learning and Computation*, 2(3), e2024JH000517.

Murray, R. J. (1996). Explicit generation of orthogonal grids for ocean models. *Journal of Computational Physics*, 126(2), 251-273.

National Energy Research Scientific Computing Center. (2025). Perlmutter architecture. Retrieved from https://docs.nersc.gov/systems/perlmutter/architecture/

Pacanowski, R. C., & Philander, S. G. H. (1981). Parameterization of vertical mixing in numerical models of tropical oceans. *Journal of Physical Oceanography*, 11(11), 1443-1451.

Paszke, A., Gross, S., Massa, F., Lerer, A., Bradbury, J., Chanan, G., ... others (2019). Pytorch: An imperative style, high-performance deep learning library. *Advances in Neural Information Processing Systems*, 32.

Paulson, C. A. (1970). The mathematical representation of wind speed and temperature profiles in the unstable atmospheric surface layer. *Journal of Applied Meteorology (1962-1982)*, 857-861.

Pawlak, G., MacCready, P., Edwards, K., & McCabe, R. (2003). Observations on the evolution of tidal vorticity at a stratified deep water headland. *Geophysical Research Letters*, 30(24), 2234. doi: 10.1029/2003GL018092

Pearson, B. (2018). Turbulence-induced anti-Stokes flow and the resulting limitations of large-eddy simulation. *Journal of Physical Oceanography*, 48(1), 117-122.

Pelletier, C., Lemarie, F., Blayo, E., Bouin, M.-N., & Redelsperger, J.-L. (2021). Two-sided turbulent surface-layer parameterizations for computing air-sea fluxes. *Quarterly Journal of the Royal Meteorological Society*, 147(736), 1726-1751.

Perez, F., & Granger, B. E. (2007). IPython: a system for interactive scientific computing. *Computing in science & engineering*, 9(3), 21-29.

Petersen, M. R., Jacobsen, D. W., Ringler, T. D., Hecht, M. W., & Maltrud, M. E. (2015). Evaluation of the arbitrary Lagrangian-Eulerian vertical coordinate method in the MPAS-Ocean model. *Ocean Modelling*, 86, 93-113.

Phillips, N. A. (1956). The general circulation of the atmosphere: A numerical experiment. *Quarterly Journal of the Royal Meteorological Society*, 82(352), 123-164.

Pressel, K. G., Kaul, C. M., Schneider, T., Tan, Z., & Mishra, S. (2015). Large-eddy simulation in an anelastic framework with closed water and entropy balances. *Journal of Advances in Modeling Earth Systems*, 7(3), 1425-1456.

Pressel, K. G., Mishra, S., Schneider, T., Kaul, C. M., & Tan, Z. (2017). Numerics and subgrid-scale modeling in large eddy simulations of stratocumulus clouds. *Journal of Advances in Modeling Earth Systems*, 9(2), 1342-1365.

Ramadhan, A., Wagner, G., Hill, C., Campin, J.-M., Churavy, V., Besard, T., ... Marshall, J. (2020). Oceananigans.jl: Fast and friendly geophysical fluid dynamics on GPUs. *Journal of Open Source Software*, 5(53). doi: 10.21105/joss.02018

Rathgeber, F., Ham, D. A., Mitchell, L., Lange, M., Luporini, F., McRae, A. T., ... Kelly, P. H. (2016). Firedrake: automating the finite element method by composing abstractions. *ACM Transactions on Mathematical Software (TOMS)*, 43(3), 1-27.

Redi, M. H. (1982). Oceanic isopycnal mixing by coordinate rotation. *Journal of Physical Oceanography*, 12(10), 1154-1158.

Reichl, B. G., & Li, Q. (2019). A parameterization with a constrained potential energy conversion rate of vertical mixing due to Langmuir turbulence. *Journal of Physical Oceanography*, 49(11), 2935-2959.

Richardson, M. I., Toigo, A. D., & Newman, C. E. (2007). PlanetWRF: A general purpose, local to global numerical model for planetary atmospheric and climate dynamics. *Journal of Geophysical Research: Planets*, 112(E9).

Ringler, T., Petersen, M., Higdon, R. L., Jacobsen, D., Jones, P. W., & Maltrud, M. (2013). A multi-resolution approach to global ocean modeling. *Ocean Modelling*, 69, 211-232.

Roquet, F., Madec, G., Brodeau, L., & Nycander, J. (2015). Defining a simplified yet "realistic" equation of state for seawater. *Journal of Physical Oceanography*, 45(10), 2564-2579.

Roquet, F., Madec, G., McDougall, T. J., & Barker, P. M. (2015). Accurate polynomial expressions for the density and specific volume of seawater using the TEOS-10 standard. *Ocean Modelling*, 90, 29-43.

Rozema, W., Bae, H. J., Moin, P., & Verstappen, R. (2015). Minimum-dissipation models for large-eddy simulation. *Physics of Fluids*, 27(8).

Schneider, T., Lan, S., Stuart, A., & Teixeira, J. (2017). Earth system modeling 2.0: A blueprint for models that learn from observations and targeted high-resolution simulations. *Geophysical Research Letters*, 44(24), 12-396.

Schumann, U., & Sweet, R. A. (1988). Fast Fourier transforms for direct solution of Poisson's equation with staggered boundary conditions. *Journal of Computational Physics*, 75(1), 123-137.

Shaw, J., & Weller, H. (2016). Comparison of terrain-following and cut-cell grids using a nonhydrostatic model. *Monthly Weather Review*, 144(6), 2085-2099.

Shchepetkin, A. F., & McWilliams, J. C. (2005). The regional oceanic modeling system (ROMS): a split-explicit, free-surface, topography-following-coordinate oceanic model. *Ocean modelling*, 9(4), 347-404.

Shu, C. (1997). Essentially non-oscillatory and weighted essentially non-oscillatory schemes for hyperbolic conservation laws (ICASE Report No. 97-65). Institute for Computer Applications in Science and Engineering, NASA Langley Research Center.

Shu, C.-W. (2009). High order weighted essentially nonoscillatory schemes for convection dominated problems. *SIAM review*, 51(1), 82-126.

Silvestri, S., Wagner, G. L., Campin, J.-M., Constantinou, N. C., Hill, C. N., Souza, A., & Ferrari, R. (2024). A new WENO-based momentum advection scheme for simulations of ocean mesoscale turbulence. *Journal of Advances in Modeling Earth Systems*, 16(7), e2023MS004130. doi: 10.1029/2023MS004130

Silvestri, S., Wagner, G. L., Constantinou, N. C., Hill, C. N., Campin, J.-M., Souza, A. N., ... Ferrari, R. (2025). A GPU-based ocean dynamical core for routine mesoscale-resolving climate simulations. *Journal of Advances in Modeling Earth Systems*, 17(4), e2024MS004465.

Silvestri, S., Wagner, G. L., Hill, C., Ardakani, M. R., Blaschke, J., Campin, J.-M., ... Ferrari, R. (2023). Oceananigans.jl: A model that achieves breakthrough resolution, memory and energy efficiency in global ocean simulations. *arXiv preprint arXiv:2309.06662*. doi: 10.48550/arXiv.2309.06662

Smagorinsky, J. (1963). General circulation experiments with the primitive equations: I. The basic experiment. *Monthly weather review*, 91(3), 99-164.

Solvik, K., Penny, S. G., & Hoyer, S. (2025). 4D-Var using Hessian approximation and backpropagation applied to automatically differentiable numerical and machine learning models. *Journal of Advances in Modeling Earth Systems*, 17(4), e2024MS004608.

Stone, P. H. (1971). Baroclinic stability under non-hydrostatic conditions. *Journal of Fluid Mechanics*, 45(4), 659-671.

Strong-Wright, J., Chen, S., Constantinou, N. C., Silvestri, S., Wagner, G. L., & Taylor, J. R. (2023). OceanBioME.jl: A flexible environment for modelling the coupled interactions between ocean biogeochemistry and physics. *Journal of Open Source Software*, 8(90), 5669.

Tsujino, H., Urakawa, S., Nakano, H., Small, R. J., Kim, W. M., Yeager, S. G., ... others (2018). JRA-55 based surface dataset for driving ocean-sea-ice models (JRA55-do). *Ocean Modelling*, 130, 79-139.

Umlauf, L., & Burchard, H. (2003). A generic length-scale equation for geophysical turbulence models. *Journal of Marine Research*, 61.

Umlauf, L., & Burchard, H. (2005). Second-order turbulence closure models for geophysical boundary layers. A review of recent work. *Continental Shelf Research*, 25(7-8), 795-827.

Vanneste, J., & Young, W. R. (2022). Stokes drift and its discontents. *Philosophical Transactions of the Royal Society A*, 380(2225), 20210032.

Vreugdenhil, C. A., & Taylor, J. R. (2018). Large-eddy simulations of stratified plane Couette flow using the anisotropic minimum-dissipation model. *Physics of Fluids*, 30(8).

Wagner, G. L., Chini, G. P., Ramadhan, A., Gallet, B., & Ferrari, R. (2021). Near-inertial waves and turbulence driven by the growth of swell. *Journal of Physical Oceanography*, 51(5), 1337-1351.

Wagner, G. L., & Constantinou, N. C. (2025). Phenomenology of decaying turbulence beneath surface waves. *Journal of Fluid Mechanics*, 1020, A51. doi: 10.1017/jfm.2025.10649

Wagner, G. L., Hillier, A., Constantinou, N. C., Silvestri, S., Souza, A., Burns, K. J., ... Ferrari, R. (2025). Formulation and calibration of CATKE, a one-equation parameterization for microscale ocean mixing. *Journal of Advances in Modeling Earth Systems*, 17(4), e2024MS004522.

Wagner, G. L., Silvestri, S., Constantinou, N. C., Strong-Wright, J., Sohail, T., Byrne, S., ... Shen, Z. (2025). CliMA/ClimaOcean.jl: v0.8.7. *Zenodo*. doi: 10.5281/zenodo.17357765

Warner, S. J., & MacCready, P. (2014). The dynamics of pressure and form drag on a sloping headland: Internal waves versus eddies. *Journal of Geophysical Research: Oceans*, 119(3), 1554-1571.

Watt-Meyer, O., Dresdner, G., McGibbon, J., Clark, S. K., Henn, B., Duncan, J., ... others (2023). ACE: A fast, skillful learned global atmospheric model for climate prediction. *arXiv preprint arXiv:2310.02074*.

Watt-Meyer, O., Henn, B., McGibbon, J., Clark, S. K., Kwa, A., Perkins, W. A., ... Bretherton, C. S. (2024). ACE2: Accurately learning subseasonal to decadal atmospheric variability and forced responses. *arXiv preprint arXiv:2411.11268*.

Yamazaki, H., Satomura, T., & Nikiforakis, N. (2016). Three-dimensional cut-cell modelling for high-resolution atmospheric simulations. *Quarterly Journal of the Royal Meteorological Society*, 142(696), 1335-1350.

Yatunin, D., Byrne, S., Kawczynski, C., Kandala, S., Bozzola, G., Sridhar, A., ... others (2025). The Climate Modeling Alliance atmosphere dynamical core: Concepts, numerics, and scaling. *Authorea Preprints*. doi: 10.22541/essoar.173940262.23304403/v1

Young, W. R. (2012). An exact thickness-weighted average formulation of the Boussinesq equations. *Journal of Physical Oceanography*, 42(5), 692-707.

Yuval, J., Langmore, I., Kochkov, D., & Hoyer, S. (2024). Neural general circulation models optimized to predict satellite-based precipitation observations. *arXiv preprint arXiv:2412.11973*.

Zhang, W., Kuo, Y.-H., Silvestri, S., Adcroft, A., Hallberg, R. W., & Griffies, S. M. (2025). A WENO finite-volume scheme for the evolution of potential vorticity in isopycnal ocean models. *Authorea Preprints*. doi: 10.22541/essoar.175380391.18723979/v1
