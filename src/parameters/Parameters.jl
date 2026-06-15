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
const AbstractEmptyParams = Union{AbstractParameters, Nothing}

"""
    mutable struct Parameters{PPHY <: AbstractEmptyParams, PSIM <: AbstractEmptyParams, PHY <: AbstractEmptyParams,
        PSOL <: AbstractEmptyParams, PUDE <: AbstractEmptyParams}

A mutable struct that holds various parameter sets for different aspects of a simulation or model.

# Fields

  - `physical::PPHY`: Physical parameters.
  - `simulation::PSIM`: Simulation parameters.
  - `hyper::PHY`: Hyperparameters.
  - `solver::PSOL`: Solver parameters.
  - `UDE::PUDE`: Universal Differential Equation (UDE) parameters.

# Type Parameters

  - `PPHY`: Type of the physical parameters, must be a subtype of `AbstractEmptyParams`.
  - `PSIM`: Type of the simulation parameters, must be a subtype of `AbstractEmptyParams`.
  - `PHY`: Type of the hyperparameters, must be a subtype of `AbstractEmptyParams`.
  - `PSOL`: Type of the solver parameters, must be a subtype of `AbstractEmptyParams`.
  - `PUDE`: Type of the UDE parameters, must be a subtype of `AbstractEmptyParams`.
"""
mutable struct Parameters{
    PPHY <: AbstractEmptyParams, PSIM <: AbstractEmptyParams, PHY <: AbstractEmptyParams,
    PSOL <: AbstractEmptyParams, PUDE <: AbstractEmptyParams}
    physical::PPHY
    simulation::PSIM
    hyper::PHY
    solver::PSOL
    UDE::PUDE
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
        simulation::SimulationParameters = SimulationParameters()
)

    # Build the parameters based on all the subtypes of parameters
    parameters = Parameters{typeof(physical), typeof(simulation), Nothing, Nothing, Nothing}(
        physical, simulation,
        nothing, nothing, nothing)

    enable_multiprocessing(parameters.simulation.multiprocessing ?
                           parameters.simulation.workers : 0)

    return parameters
end

function Base.:(==)(a::Parameters, b::Parameters)
    a.physical == b.physical && a.simulation == b.simulation &&
        a.solver == b.solver && a.hyper == b.hyper &&
        a.UDE == b.UDE
end

# Display setup
Base.show(io::IO, ::MIME"text/plain", params::Parameters) = Base.show(io, params)
function Base.show(io::IO, params::Parameters)
    pad = 14

    println(io, "Parameters")

    # Physical
    label(io, "  physical", pad)
    if isnothing(params.physical)
        hint(io, "(nothing)")
    else
        p = params.physical
        field(io, "ρ");
        print(io, " = ");
        val(io, "$(p.ρ)")
        sep(io)
        field(io, "A");
        print(io, " ∈ [");
        val(io, "$(p.minA)");
        print(io, ", ");
        val(io, "$(p.maxA)");
        print(io, "]")
        sep(io)
        field(io, "C");
        print(io, " ∈ [");
        val(io, "$(p.minC)");
        print(io, ", ");
        val(io, "$(p.maxC)");
        print(io, "]")
    end
    println(io)

    # Simulation
    label(io, "  simulation", pad)
    if isnothing(params.simulation)
        hint(io, "(nothing)")
    else
        s = params.simulation
        field(io, "tspan");
        print(io, " = ");
        val(io, "$(s.tspan)")
        sep(io)
        print(io, check(s.use_iceflow));
        field(io, "iceflow")
        print(io, " ")
        print(io, check(s.use_MB));
        field(io, "MB")
        print(io, " ")
        print(io, check(s.use_velocities));
        field(io, "velocities")
    end
    println(io)

    # Solver
    label(io, "  solver", pad)
    if isnothing(params.solver)
        hint(io, "(nothing)")
    else
        sv = params.solver
        val(io, "$(nameof(typeof(sv.solver)))")
        sep(io)
        field(io, "reltol");
        print(io, " = ");
        val(io, "$(sv.reltol)")
        sep(io)
        field(io, "maxiters");
        print(io, " = ");
        val(io, "$(sv.maxiters)")
    end
    println(io)

    # Hyper
    label(io, "  hyper", pad)
    if isnothing(params.hyper)
        hint(io, "(nothing)")
    else
        h = params.hyper
        field(io, "epochs");
        print(io, " = ");
        val(io, "$(h.epochs)")
        sep(io)
        field(io, "batch_size");
        print(io, " = ");
        val(io, "$(h.batch_size)")
        sep(io)
        field(io, "optimizer");
        print(io, " = ")
        if h.optimizer isa Vector
            val(io, "[")
            first = true
            for o in h.optimizer
                first || val(io, ", ")
                first = false
                val(io, nameof(typeof(o)))
            end
            val(io, "]")
        else
            val(io, "$(nameof(typeof(h.optimizer)))")
        end
    end
    println(io)

    # UDE
    label(io, "  UDE", pad)
    if isnothing(params.UDE)
        hint(io, "(nothing)")
    else
        u = params.UDE
        field(io, "target");
        print(io, " = ")
        isnothing(u.target) ? hint(io, "(nothing)") : val(io, ":$(u.target)")
        sep(io)
        val(io, "\"$(u.optimization_method)\"")
        sep(io)
        isnothing(u.grad) ? hint(io, "(nothing)") : val(io, "$(nameof(typeof(u.grad)))")
        sep(io)
        field(io, "loss");
        print(io, " = ");
        val(io, "$(nameof(typeof(u.empirical_loss_function)))")
    end
    println(io)
end
