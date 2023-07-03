import Pkg
Pkg.activate(dirname(Base.current_project()))

using Revise
using Sleipnir
using Test
using JLD2
using Infiltrator

include(joinpath(Sleipnir.root_dir, "test/params_construction.jl"))

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

@testset "Parameters constructors with specified values" params_constructor_specified(save_refs=false)

@testset "Parameters constructors by default" params_constructor_default(save_refs=false)

