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
using Rasters

# ##############################################
# ############    PARAMETERS     ###############
# ##############################################

cd(@__DIR__)
const global root_dir::String = dirname(Base.current_project())

# ##############################################
# ##########  SLEIPNIR LIBRARIES  ##############
# ##############################################

include("setup/config.jl")
# All parameters needed for the models
include("parameters/Parameters.jl")
# Anything related to managing glacier topographical and climate data
include("glaciers/glacier/Glacier.jl")
# All structures and functions related to ODINN models
include("models/Model.jl")
# Everything related to running simulations in ODINN
include("simulations/Simulation.jl")

end # module
