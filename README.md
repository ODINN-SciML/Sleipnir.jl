[![Build Status](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ODINN-SciML/Sleipnir.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ODINN-SciML/Sleipnir.jl)
[![CompatHelper](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CompatHelper.yml)
[![docs](https://img.shields.io/badge/documentation-dev-blue.svg)](https://odinn-sciml.github.io/ODINN.jl/dev/Packages/sleipnir/)

<img src="https://github.com/JordiBolibar/Sleipnir.jl/blob/main/data/Sleipnir_logo-19.png" width="250">

## About Sleipnir.jl

Sleipnir.jl is the core package of [ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl), containing all the basic data structures to manage glacier and climate data, as well as multiple types of numerical simulations and parameters.

It provides:

  - Glacier and climate containers: `Glacier2D`, `Climate2D`, and the observation data attached to them (ice thickness, surface velocities, geodetic elevation change).
  - The parameter hierarchy: `Parameters`, `SimulationParameters`, `PhysicalParameters`.
  - The law abstraction (`Law`, `AbstractLaw`) used to plug physical or machine-learning computations into the ice flow and mass balance models, including the VJP infrastructure needed for inverse modelling.
  - The `Model` container and the `Results` container with post-processing and plotting utilities.

Sleipnir is part of the ODINN ecosystem, where each package has a narrow role:

  - [Gungnir](https://github.com/ODINN-SciML/Gungnir) (Python): preprocesses OGGM glacier and climate data, read by Sleipnir.
  - **Sleipnir**: core data structures (this package). All the other Julia packages depend on it.
  - [Muninn.jl](https://github.com/ODINN-SciML/Muninn.jl): surface mass balance models.
  - [Huginn.jl](https://github.com/ODINN-SciML/Huginn.jl): ice flow models and PDE solvers.
  - [ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl): differentiable pipeline for UDE training and inversions.

## Use Sleipnir directly or ODINN.jl?

Most users should install [ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl), which re-exports everything in Sleipnir. Use Sleipnir on its own when you want to build or inspect glacier and climate data structures without loading the whole simulation stack, or when you are prototyping a new `Law` or dynamic input that will later be used in Huginn or ODINN.

## Installing Sleipnir

> `Sleipnir.jl` requires Julia v1.11.

In order to install `Sleipnir` in a given environment, just do in the REPL:
```julia
julia> ] # enter Pkg mode
(@v1.11) pkg> activate MyEnvironment # or activate whatever path for the Julia environment
(MyEnvironment) pkg> add Sleipnir
```

The preprocessed glacier data are downloaded automatically the first time Sleipnir is precompiled (see [Data preprocessing](#data-preprocessing)).

## How to use Sleipnir

The following example loads one glacier and inspects its initial state:

```julia
using Sleipnir

# Multiprocessing is disabled for local runs
params = Parameters(
    simulation = SimulationParameters(
        tspan = (2010.0, 2015.0),
        multiprocessing = false,
        rgi_paths = get_rgi_paths()
    )
)

# Initialize the glacier from the preprocessed data
glaciers = initialize_glaciers(["RGI60-11.03638"], params)
glacier = glaciers[1]

@show glacier.rgi_id, glacier.nx, glacier.ny
@show size(glacier.H₀) # initial ice thickness
@show size(glacier.S)  # surface elevation
```

To run simulations on these glaciers, see the tutorials in the [ODINN documentation](https://odinn-sciml.github.io/ODINN.jl/dev/forward_simulation/). Sleipnir's own page is [here](https://odinn-sciml.github.io/ODINN.jl/dev/Packages/sleipnir/), the full list of types and functions is in the [API reference](https://odinn-sciml.github.io/ODINN.jl/dev/API/api_sleipnir/), and guidance on adding new laws, inputs or data is in [Extending ODINN](https://odinn-sciml.github.io/ODINN.jl/dev/extending/).

## Data preprocessing

As of version 0.7.1, OGGM data are now preprocessed with [Gungnir](https://github.com/ODINN-SciML/Gungnir). These preprocessed data are saved on a Hugging Face repository they are downloaded as artifacts upon precompilation of Sleipnir. They are then stored locally in `~/.ODINN/ODINN_prepro/` for the subsequent executions.

In case for example you want to perform simulations with glaciers that are not in the preprocessed directory, the preprocessed directory path can be overridden very easily.
To do this, define an `Overrides.toml`, which should be placed in `~/.julia/artifacts/Overrides.toml`.
It must contain the UUID of Sleipnir together with the path to your custom preprocessed directory:
```
[f5e6c550-199f-11ee-3608-394420200519]
ODINN_prepro = "/path/to/custom/dir"
```
See [the artifacts documentation](https://pkgdocs.julialang.org/v1/artifacts/) for more information.

## Contributing and community

Contributions are welcome. You can report bugs and request features in the [issues](https://github.com/ODINN-SciML/Sleipnir.jl/issues) tab, or open a pull request against `main` from a fork. See [How to contribute](https://odinn-sciml.github.io/ODINN.jl/dev/contribute/) and the [Code of conduct](https://odinn-sciml.github.io/ODINN.jl/dev/code_of_conduct/) for the guidelines shared across the ODINN ecosystem.

## How to cite

If you use Sleipnir, please cite the ODINN paper published in [Geoscientific Model Development](https://gmd.copernicus.org/articles/16/6671/2023/gmd-16-6671-2023.html):
```
@article{bolibar_sapienza_universal_2023,
	title = {Universal differential equations for glacier ice flow modelling},
	author = {Bolibar, J. and Sapienza, F. and Maussion, F. and Lguensat, R. and Wouters, B. and P\'erez, F.},
	journal = {Geoscientific Model Development},
	volume = {16},
	year = {2023},
	number = {22},
	pages = {6671--6687},
	url = {https://gmd.copernicus.org/articles/16/6671/2023/},
	doi = {10.5194/gmd-16-6671-2023}
}
```
