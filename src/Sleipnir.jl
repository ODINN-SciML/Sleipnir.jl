__precompile__() # this module is safe to precompile
module Sleipnir

# ##############################################
# ###########       PACKAGES     ##############
# ##############################################

using Base: @kwdef
using Infiltrator
import Pkg
using PyCall

# ##############################################
# ############    PARAMETERS     ###############
# ##############################################

const global root_dir::String = dirname(Base.current_project())


# ##############################################
# ############  PYTHON LIBRARIES  ##############
# ##############################################

# const netCDF4::PyObject = PyNULL()
# const cfg::PyObject = PyNULL()
# const utils::PyObject = PyNULL()
# const workflow::PyObject = PyNULL()
# const tasks::PyObject = PyNULL()
# const global_tasks::PyObject = PyNULL()
# const graphics::PyObject = PyNULL()
# const bedtopo::PyObject = PyNULL()
# const millan22::PyObject = PyNULL()
# const MBsandbox::PyObject = PyNULL()
# const salem::PyObject = PyNULL()

# # Essential Python libraries
# const np::PyObject = PyNULL()
# const xr::PyObject = PyNULL()
# const rioxarray::PyObject = PyNULL()
# const pd::PyObject = PyNULL()

# ##############################################
# ##########  SLEIPNIR LIBRARIES  ##############
# ##############################################

# All parameters needed for the models
include(joinpath(Sleipnir.root_dir, "src/parameters/Parameters.jl"))
# Anything related to managing glacier topographical and climate data
include(joinpath(Sleipnir.root_dir, "src/glaciers/Glacier.jl"))
# All structures and functions related to ODINN models
include(joinpath(Sleipnir.root_dir, "src/models/Model.jl"))
# Everything related to running simulations in ODINN
include(joinpath(Sleipnir.root_dir, "src/simulations/Simulation.jl"))

end # module