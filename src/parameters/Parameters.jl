export Parameters, PhysicalParameters, SimulationParameters

abstract type AbstractParameters end

include("PhysicalParameters.jl")
include("SimulationParameters.jl")

struct Parameters{P <: AbstractParameters}
    simulation::Union{P,Nothing}
    physical::Union{P,Nothing}
    hyper::Union{P,Nothing}
    solver::Union{P,Nothing}
    UDE::Union{P,Nothing} 
    OGGM::Union{P,Nothing}
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
            simulation::SimulationParameters = SimulationParameters(),
            physical::PhysicalParameters = PhysicalParameters()
            )

    # Build the parameters based on all the subtypes of parameters
    parameters = Parameters(simulation, physical, nothing, nothing,
                            nothing, nothing)

    return parameters
end
