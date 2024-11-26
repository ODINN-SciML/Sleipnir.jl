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

# ##############################################
# ############    PARAMETERS     ###############
# ##############################################

cd(@__DIR__)
const global root_dir::String = dirname(Base.current_project())

# ##############################################
# ############  PYTHON LIBRARIES  ##############
# ##############################################

# We first load the SSL package 
# @eval using OpenSSL_jll
using PythonCall, CondaPkg
# @eval using PythonCall, CondaPkg
# const openssl = Ref{Py}()
# @eval const openssl = pyimport("ssl")

# We define empty objects for the Python packages
const netCDF4 = Ref{Py}()
const cfg = Ref{Py}()
const utils = Ref{Py}()
const workflow = Ref{Py}()
const tasks = Ref{Py}()
const global_tasks = Ref{Py}()
const graphics = Ref{Py}()
const bedtopo = Ref{Py}()
const millan22 = Ref{Py}()
const MBsandbox = Ref{Py}()
const salem = Ref{Py}()

# Essential Python libraries
const xr = Ref{Py}()
const rioxarray = Ref{Py}()
const pd = Ref{Py}()

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
