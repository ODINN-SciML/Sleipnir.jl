import Pkg
Pkg.activate(dirname(Base.current_project()))

if !parse(Bool, get(ENV, "CI", "false"))
    using Revise
end
using Sleipnir
using Test
using JLD2
using Infiltrator
using CairoMakie
using JSON
using Rasters
import NCDatasets
using DimensionalData
using Dates
using DateFormats
using JET

include("params_construction.jl")
include("glaciers_construction.jl")
include("surface_velocity.jl")
include("plot_utils.jl")
include("results.jl")

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

@testset "Run all tests" begin

@testset "Constructors" begin
    @testset "Parameters constructors with specified values" params_constructor_specified(save_refs=false)
    @testset "Parameters constructors by default" params_constructor_default(save_refs=false)
    @testset "Glaciers 2D constructors w/o glathida data" glaciers2D_constructor(use_glathida_data=false, save_refs=false)
    @testset "Glaciers 2D constructors w/ glathida data" glaciers2D_constructor(use_glathida_data=true, save_refs=false)
    @testset "Surface velocity datacube" surface_velocity_data()
    @testset "Results instantiation w/o velocity datacube" results_default(save_refs=false)
    @testset "Results instantiation w/ velocity datacube" results_default(save_refs=false, useDatacube=true)
end

@testset "Plotting functions" begin
    @testset "Glaciers 2D plots" glaciers2D_plots()
    @testset "Video plot test" make_thickness_video_test()
end

end
