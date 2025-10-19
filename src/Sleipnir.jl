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
using Observables
import Contour
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
using GR
using DataStructures
using ImageInTerminal
using ImageCore
using Printf
using CFTime
using MLStyle
using FFTW
using Zygote # To skip some lines in Sleipnir/src/simulations/results/results_utils.jl that causes an error with SciMLSensitivity

##############################################
############    PARAMETERS     ###############
##############################################

cd(@__DIR__)
const global root_dir::String = dirname(Base.current_project())
const global prepro_dir::String = joinpath(homedir(), ".ODINN", "ODINN_prepro")
const doublePrec::Bool = parse(Bool, get(ENV, "ODINN_DOUBLE_PREC", "true"))
const Float = doublePrec ? Float64 : Float32
const Int = doublePrec ? Int64 : Int32
if !doublePrec
    @warn "Double precision is disabled"
end

##############################################
##########  SLEIPNIR LIBRARIES  ##############
##############################################

include("setup/config.jl")

# Anything related to managing glacier data used for data assimilation
include("glaciers/data/Data.jl")

# All parameters needed for the models
include("parameters/Parameters.jl")

# Anything related to managing glacier topographical and climate variables
include("glaciers/glacier/Glacier.jl")

# The utils of surface velocity data, glaciers and climate need the struct to be already
# defined since they depend on each other. This is why we import them afterwards
include("glaciers/data/SurfaceVelocityData_utils.jl")
include("glaciers/data/SurfaceVelocityMapping_utils.jl")
include("glaciers/glacier/glacier2D_utils.jl")
include("glaciers/climate/climate2D_utils.jl")

# All structures and functions related to ODINN models
include("models/Model.jl")

# Everything related to running simulations in ODINN
include("simulations/Simulation.jl")
# Law interface and utils
include("laws/Cache.jl")
include("laws/Law.jl")

# Fake data used in the tests
include("data/surface_velocity.jl")

##############################################
#######    PRE-LOADED VARIABLES     ##########
##############################################

end # module
