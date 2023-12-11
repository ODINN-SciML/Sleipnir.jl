__precompile__() # this module is safe to precompile
module Sleipnir

# ##############################################
# ###########       PACKAGES     ##############
# ##############################################

using Base: @kwdef
using Infiltrator
import Pkg
using PyCall
using JLD2
using Distributed
using Statistics
using CairoMakie
using Downloads
using HDF5

# ##############################################
# ############    PARAMETERS     ###############
# ##############################################

cd(@__DIR__)
const global root_dir::String = dirname(Base.current_project())

# ##############################################
# ############  PYTHON LIBRARIES  ##############
# ##############################################

# We either retrieve the reexported Python libraries from Sleipnir or we start from scratch
const netCDF4::PyObject = isdefined(Sleipnir, :netCDF4) ? Sleipnir.netCDF4 : PyNULL()
const cfg::PyObject = isdefined(Sleipnir, :cfg) ? Sleipnir.cfg : PyNULL()
const utils::PyObject = isdefined(Sleipnir, :utils) ? Sleipnir.utils : PyNULL()
const workflow::PyObject = isdefined(Sleipnir, :workflow) ? Sleipnir.workflow : PyNULL()
const tasks::PyObject = isdefined(Sleipnir, :tasks) ? Sleipnir.tasks : PyNULL()
const global_tasks::PyObject = isdefined(Sleipnir, :global_tasks) ? Sleipnir.global_tasks : PyNULL()
const graphics::PyObject = isdefined(Sleipnir, :graphics) ? Sleipnir.graphics : PyNULL()
const bedtopo::PyObject = isdefined(Sleipnir, :bedtopo) ? Sleipnir.bedtopo : PyNULL()
const millan22::PyObject = isdefined(Sleipnir, :millan22) ? Sleipnir.millan22 : PyNULL()
const MBsandbox::PyObject = isdefined(Sleipnir, :MBsandbox) ? Sleipnir.MBsandbox : PyNULL()
const salem::PyObject = isdefined(Sleipnir, :salem) ? Sleipnir.salem : PyNULL()

# Essential Python libraries
const xr::PyObject = isdefined(Sleipnir, :xr) ? Sleipnir.xr : PyNULL()
const rioxarray::PyObject = isdefined(Sleipnir, :rioxarray) ? Sleipnir.rioxarray : PyNULL()
const pd::PyObject = isdefined(Sleipnir, :pd) ? Sleipnir.pd : PyNULL()

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
