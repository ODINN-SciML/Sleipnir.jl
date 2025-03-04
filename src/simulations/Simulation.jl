
export Simulation, Results

"""
    Simulation

An abstract type representing a generic simulation. This type is intended to be 
subclassed by specific simulation types to provide a common interface and shared 
functionality for all simulations.
"""
abstract type Simulation end

include("results/Results.jl")

###############################################
################### UTILS #####################
###############################################

include("simulation_utils.jl")
include("results/results_utils.jl")
include("results/results_plotting_utils.jl")
include("results/results_plotting_video_utils.jl")
