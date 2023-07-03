
# Abstract type as a parent type for simulations
abstract type Simulation end

include("results/Results.jl")

###############################################
################### UTILS #####################
###############################################

include(joinpath(Sleipnir.root_dir, "src/simulations/simulation_utils.jl"))
include(joinpath(Sleipnir.root_dir, "src/simulations/results/results_utils.jl"))


