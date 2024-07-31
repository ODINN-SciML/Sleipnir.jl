import Pkg
Pkg.activate(dirname(Base.current_project()))

using Revise
using Sleipnir
using Test
using JLD2
using Infiltrator
using CairoMakie

include("params_construction.jl")
include("glaciers_construction.jl")
include("plot_utils.jl")

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

@testset "Parameters constructors with specified values" params_constructor_specified()

@testset "Parameters constructors by default" params_constructor_default()

@testset "Glaciers 2D constructors" glaciers2D_constructor()

#@testset "Glaciers 2D plots" glaciers2D_plots()