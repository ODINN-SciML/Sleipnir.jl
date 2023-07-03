export Parameters, AbstractParameters, PhysicalParameters, SimulationParameters

abstract type AbstractParameters end

const AbstractEmptyParams = Union{AbstractParameters,Nothing}

include("PhysicalParameters.jl")
include("SimulationParameters.jl")

struct Parameters{PPHY <: AbstractEmptyParams, PSIM <: AbstractEmptyParams, PHY <: AbstractEmptyParams, 
                  PSOL <: AbstractEmptyParams, PUDE <: AbstractEmptyParams, POGGM <: AbstractEmptyParams}
    physical::PPHY
    simulation::PSIM
    hyper::PHY
    solver::PSOL
    UDE::PUDE
    OGGM::POGGM
end

"""
Parameters(;
        simulation::SimulationParameters = SimulationParameters()
        physical::PhysicalParameters = PhysicalParameters()
        )
Initialize ODINN parameters

Keyword arguments
=================
    
"""
function Parameters(;
            physical::PhysicalParameters = PhysicalParameters(),
            simulation::SimulationParameters = SimulationParameters()
            ) 

    # Build the parameters based on all the subtypes of parameters
    parameters = Parameters(physical, simulation, 
                            nothing, nothing,nothing, nothing)

    return parameters
end
