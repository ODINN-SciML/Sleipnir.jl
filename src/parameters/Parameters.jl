export AbstractParameters, PhysicalParameters, SimulationParameters, AbstractEmptyParams

"""
    AbstractParameters

An abstract type that serves as a base for all parameter-related types in the ODINN ecosystem.
"""
abstract type AbstractParameters end

"""
        AbstractEmptyParams

A type alias that represents a union of `AbstractParameters` and `Nothing`.
This can be used to indicate that a parameter can either be an instance of
`AbstractParameters` or `nothing`.
"""
const AbstractEmptyParams = Union{AbstractParameters,Nothing}

"""
    mutable struct Parameters{PPHY <: AbstractEmptyParams, PSIM <: AbstractEmptyParams, PHY <: AbstractEmptyParams,
        PSOL <: AbstractEmptyParams, PUDE <: AbstractEmptyParams, PINV <: AbstractEmptyParams}

A mutable struct that holds various parameter sets for different aspects of a simulation or model.

# Fields
- `physical::PPHY`: Physical parameters.
- `simulation::PSIM`: Simulation parameters.
- `hyper::PHY`: Hyperparameters.
- `solver::PSOL`: Solver parameters.
- `UDE::PUDE`: Universal Differential Equation (UDE) parameters.
- `inversion::PINV`: Inversion parameters.

# Type Parameters
- `PPHY`: Type of the physical parameters, must be a subtype of `AbstractEmptyParams`.
- `PSIM`: Type of the simulation parameters, must be a subtype of `AbstractEmptyParams`.
- `PHY`: Type of the hyperparameters, must be a subtype of `AbstractEmptyParams`.
- `PSOL`: Type of the solver parameters, must be a subtype of `AbstractEmptyParams`.
- `PUDE`: Type of the UDE parameters, must be a subtype of `AbstractEmptyParams`.
- `PINV`: Type of the inversion parameters, must be a subtype of `AbstractEmptyParams`.
"""
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
    Parameters(; physical::PhysicalParameters = PhysicalParameters(), simulation::SimulationParameters = SimulationParameters())

Constructs a `Parameters` object with the given physical and simulation parameters.

# Arguments
- `physical::PhysicalParameters`: An instance of `PhysicalParameters` (default: `PhysicalParameters()`).
- `simulation::SimulationParameters`: An instance of `SimulationParameters` (default: `SimulationParameters()`).

# Returns
- A `Parameters` object initialized with the provided physical and simulation parameters.

# Notes
- If `simulation.multiprocessing` is enabled, multiprocessing is configured with the specified number of workers.
"""
function Parameters(;
    physical::PhysicalParameters = PhysicalParameters(),
    simulation::SimulationParameters = SimulationParameters(),
)

    # Build the parameters based on all the subtypes of parameters
    parameters = Parameters(
        physical, simulation,
        nothing, nothing, nothing, nothing)

    enable_multiprocessing(parameters.simulation.multiprocessing ? parameters.simulation.workers : 0)

    return parameters
end

Base.:(==)(a::Parameters, b::Parameters) = a.physical == b.physical && a.simulation == b.simulation &&
                                           a.solver == b.solver && a.hyper == b.hyper &&
                                           a.UDE == b.UDE && a.inversion == b.inversion
