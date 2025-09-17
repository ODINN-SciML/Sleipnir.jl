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

include("iceflow_def.jl")
include("params_construction.jl")
include("glaciers_construction.jl")
include("thickness_construction.jl")
include("surface_velocity.jl")
include("climate.jl")
include("plot_utils.jl")
include("results.jl")
include("laws.jl")
include("misc.jl")

# Activate to avoid GKS backend Plot issues in the JupyterHub
ENV["GKSwstype"]="nul"

@testset "Run all tests" begin

@testset "Constructors" begin
    @testset "Parameters constructors with specified values" params_constructor_specified(save_refs=false)
    @testset "Parameters constructors by default" params_constructor_default(save_refs=false)
    @testset "Glaciers 2D constructors w/o glathida data" glaciers2D_constructor(use_glathida_data=false, save_refs=false)
    @testset "Glaciers 2D constructors w/ glathida data" glaciers2D_constructor(use_glathida_data=true, save_refs=false)
    @testset "Surface velocity datacube" surface_velocity_data()
    @testset "Thickness data constructor" thickness_construction()
    @testset "Results instantiation w/o velocity datacube" results_default(save_refs=false)
    @testset "Results instantiation w/ velocity datacube" results_default(save_refs=false, useDatacube=true)
end

@testset "Climate operations" begin
    @testset "Dummy climate" dummy_climate()
    @testset "Climate downscale" climate_downscale()
end

@testset "Misc" begin
    @testset "Helpers" helpers()
    @testset "Simulation utils" simulation_utils()
    @testset "Results utils" results_utils()
    @testset "Glacier grid downscaling" glacier_grid_downscaling()
    @testset "Operations on glacier mask" operations_glacier_mask()
end

@testset "Plotting functions" begin
    @testset "Glaciers 2D plots" glaciers2D_plots()
    @testset "Video plot test" make_thickness_video_test()
end

@testset "Law" begin
    generate_inputs_testset()
    normalize_law_inputs_testset()
    apply_law_testset()
end

end
