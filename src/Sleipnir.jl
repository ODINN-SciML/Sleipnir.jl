__precompile__() # this module is safe to precompile
module Sleipnir

# ##############################################
# ###########       PACKAGES     ##############
# ##############################################

using Base: @kwdef
using Infiltrator
import Pkg
using JLD2
using Distributed
using Statistics
using CairoMakie
using Downloads
using HDF5
using ComponentArrays
using Rasters
using CSV
using JSON
using CodecZlib
using Tar
import NCDatasets
using Unitful: m, rad, °
using CoordRefSystems
using Dates, DateFormats
using NetCDF
using GR

# ##############################################
# ############    PARAMETERS     ###############
# ##############################################

cd(@__DIR__)
const global root_dir::String = dirname(Base.current_project())
const global prepro_dir::String = joinpath(homedir(), ".ODINN", "ODINN_prepro")

# ##############################################
# ##########  SLEIPNIR LIBRARIES  ##############
# ##############################################

include("setup/config.jl")
# All parameters needed for the models
include("parameters/Parameters.jl")
# Anything related to managing glacier topographical and climate variables
include("glaciers/glacier/Glacier.jl")
# Anything related to managing glacier data used for data assimilation
include("glaciers/data/Data.jl")
# All structures and functions related to ODINN models
include("models/Model.jl")
# Everything related to running simulations in ODINN
include("simulations/Simulation.jl")

end # module