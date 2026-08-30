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
using Statistics, NaNStatistics
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
using DataStructures
using ImageInTerminal
using ImageCore
using Printf
using CFTime
using MLStyle
using FFTW
using Missings
using Zygote # To skip some lines in Sleipnir/src/simulations/results/results_utils.jl that causes an error with SciMLSensitivity

# Local meshgrid replacement (formerly from GR)
meshgrid(x, y) = (repeat(x, 1, length(y)), repeat(y', length(x), 1))

##############################################
############    PARAMETERS     ###############
##############################################

const global root_dir::String = dirname(@__DIR__)
const dir_src::String = joinpath(root_dir, "src")
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

include(dir_src*"/setup/config.jl")

# Anything related to managing glacier data used for data assimilation
include(dir_src*"/glaciers/data/Data.jl")

# All parameters needed for the models
include(dir_src*"/parameters/Parameters.jl")

# Anything related to managing glacier topographical and climate variables
include(dir_src*"/glaciers/glacier/Glacier.jl")

# The utils of surface velocity data, glaciers and climate need the struct to be already
# defined since they depend on each other. This is why we import them afterwards
include(dir_src*"/glaciers/data/SurfaceVelocityData_utils.jl")
include(dir_src*"/glaciers/data/SurfaceVelocityMapping_utils.jl")
include(dir_src*"/glaciers/glacier/glacier2D_topography.jl")
include(dir_src*"/glaciers/glacier/glacier2D_projection.jl")
include(dir_src*"/glaciers/glacier/glacier2D_geodata.jl")
include(dir_src*"/glaciers/glacier/glacier2D_utils.jl")
include(dir_src*"/glaciers/climate/climate2D_utils.jl")

# All structures and functions related to ODINN models
include(dir_src*"/models/Model.jl")

# Everything related to running simulations in ODINN
include(dir_src*"/simulations/Simulation.jl")
# Law interface and utils
include(dir_src*"/laws/GenInput.jl")
include(dir_src*"/laws/Inputs.jl")
include(dir_src*"/laws/Cache.jl")
include(dir_src*"/laws/AbstractLaw.jl")
include(dir_src*"/laws/VJP.jl")
include(dir_src*"/laws/Law.jl")

# Fake data used in the tests
include(dir_src*"/data/surface_velocity.jl")

# Abstract loss definition
include(dir_src*"/losses/Losses.jl")

##############################################
#######    PRE-LOADED VARIABLES     ##########
##############################################

end # module
