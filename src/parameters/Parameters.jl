export Parameters, AbstractParameters, PhysicalParameters, SimulationParameters, AbstractEmptyParams

abstract type AbstractParameters end

const AbstractEmptyParams = Union{AbstractParameters,Nothing}

mutable struct Parameters{PPHY <: AbstractEmptyParams, PSIM <: AbstractEmptyParams, PHY <: AbstractEmptyParams,
        PSOL <: AbstractEmptyParams, PUDE <: AbstractEmptyParams, PINV <: AbstractEmptyParams}
        physical::PPHY
        simulation::PSIM
        hyper::PHY
        solver::PSOL
        UDE::PUDE
        inversion::PINV
end

include("PhysicalParameters.jl")
include("SimulationParameters.jl")

"""
Parameters(;
        physical::PhysicalParameters = PhysicalParameters(),
        simulation::SimulationParameters = SimulationParameters()
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
                            nothing, nothing, nothing, nothing)

    if parameters.simulation.multiprocessing
            enable_multiprocessing(parameters.simulation.workers)
    end
    
    return parameters
end

Base.:(==)(a::Parameters, b::Parameters) = a.physical == b.physical && a.simulation == b.simulation && 
                                           a.solver == b.solver && a.hyper == b.hyper && 
                                           a.UDE == b.UDE && a.inversion == b.inversion  
