export Parameters, AbstractParameters, PhysicalParameters, SimulationParameters, OGGMparameters, AbstractEmptyParams

abstract type AbstractParameters end

const AbstractEmptyParams = Union{AbstractParameters,Nothing}

include("PhysicalParameters.jl")
include("SimulationParameters.jl")
include("OGGMparameters.jl")

struct Parameters{PPHY <: AbstractEmptyParams, PSIM <: AbstractEmptyParams, PHY <: AbstractEmptyParams, 
        PSOL <: AbstractEmptyParams, PUDE <: AbstractEmptyParams, POGGM <: AbstractEmptyParams, PINV <: AbstractEmptyParams}  
        physical::PPHY
        simulation::PSIM
        OGGM::POGGM
        hyper::PHY
        solver::PSOL
        UDE::PUDE
        inversion::PINV  
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
        OGGM::OGGMparameters = OGGMparameters(),
    ) 

    # Build the parameters based on all the subtypes of parameters
    parameters = Parameters(physical, simulation, OGGM,
                            nothing, nothing, nothing, nothing)  

    if parameters.simulation.multiprocessing
            enable_multiprocessing(parameters.simulation.workers)
    end
            
    oggm_config(parameters.OGGM.working_dir; oggm_processes=parameters.OGGM.workers)
    
    return parameters
end

Base.:(==)(a::Parameters, b::Parameters) = a.physical == b.physical && a.simulation == b.simulation && 
                                           a.OGGM == b.OGGM && a.solver == b.solver && a.hyper == b.hyper && 
                                           a.UDE == b.UDE && a.inversion == b.inversion  
