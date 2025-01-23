[![Build Status](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ODINN-SciML/Sleipnir.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ODINN-SciML/Sleipnir.jl)
[![CompatHelper](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CompatHelper.yml)

<img src="https://github.com/JordiBolibar/Sleipnir.jl/blob/main/data/Sleipnir_logo-19.png" width="250">

Sleipnir.jl is the core package of [ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl), containing all the basic data structures to manage glacier and climate data, as well as multiple types of numerical simulations and parameters.

## Data preprocessing

As of version 0.7.1, OGGM data are now preprocessed with [Gungnir](https://github.com/ODINN-SciML/Gungnir). These preprocessed data are saved on a server and Sleipnir downloads them during the first execution on a machine. They are then stored locally in `~/.ODINN/ODINN_prepro/` for the subsequent executions.
