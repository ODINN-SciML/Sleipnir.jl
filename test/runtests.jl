import Pkg
Pkg.activate(dirname(Base.current_project()))

using Revise
using Sleipnir
using PyCall
using Test
using JLD2
using Infiltrator

include("params_construction.jl")
include("glaciers_construction.jl")

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

@testset "Parameters constructors with specified values" params_constructor_specified()

@testset "Parameters constructors by default" params_constructor_default()

@testset "Glaciers 2D constructors" glaciers2D_constructor()
