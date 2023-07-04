export Parameters, AbstractParameters, PhysicalParameters, SimulationParameters, OGGMparameters

abstract type AbstractParameters end

const AbstractEmptyParams = Union{AbstractParameters,Nothing}

include("PhysicalParameters.jl")
include("SimulationParameters.jl")
include("OGGMparameters.jl")

struct Parameters{PPHY <: AbstractEmptyParams, PSIM <: AbstractEmptyParams, PHY <: AbstractEmptyParams, 
                  PSOL <: AbstractEmptyParams, PUDE <: AbstractEmptyParams, POGGM <: AbstractEmptyParams}
    physical::PPHY
    simulation::PSIM
    OGGM::POGGM
    hyper::PHY
    solver::PSOL
    UDE::PUDE
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
            simulation::SimulationParameters = SimulationParameters(),
            OGGM::OGGMparameters = OGGMparameters()
            ) 

    # Build the parameters based on all the subtypes of parameters
    parameters = Parameters(physical, simulation, OGGM,
                            nothing,nothing, nothing)

    return parameters
end
