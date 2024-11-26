[![Build Status](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ODINN-SciML/Sleipnir.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/ODINN-SciML/Sleipnir.jl)
[![CompatHelper](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CompatHelper.yml/badge.svg)](https://github.com/ODINN-SciML/Sleipnir.jl/actions/workflows/CompatHelper.yml)

<img src="https://github.com/JordiBolibar/Sleipnir.jl/blob/main/data/Sleipnir_logo-19.png" width="250">

Sleipnir.jl is the core package of [ODINN.jl](https://github.com/ODINN-SciML/ODINN.jl), containing all the basic data structures to manage glacier and climate data, as well as multiple types of numerical simulations and parameters.

## Installation

Sleipnir.jl (and all the libraries in ODINN.jl) use Python dependencies that are handled by CondaPkg. 
In order to these dependencies to work in harmony with the Julia package dependencies, we need to set the SSL version inside the Julia session. 
Before importing Sleipnir, we need to set the version of OpenSSL with the following command. 
```julia
using PreferenceTools
pkg > preference add CondaPkg openssl_version=3.4.0
```
And then we can simply import Sleipnir 
with 
```
using Sleipnir
```
The first time you compile Sleipir, this will install the new conda package in `.CondaPkg/`.