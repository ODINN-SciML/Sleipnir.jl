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
mutable struct Parameters{
    PPHY <: AbstractEmptyParams, PSIM <: AbstractEmptyParams, PHY <: AbstractEmptyParams,
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
        simulation::SimulationParameters = SimulationParameters()
)

    # Build the parameters based on all the subtypes of parameters
    parameters = Parameters(
        physical, simulation,
        nothing, nothing, nothing, nothing)

    enable_multiprocessing(parameters.simulation.multiprocessing ?
                           parameters.simulation.workers : 0)

    return parameters
end

function Base.:(==)(a::Parameters, b::Parameters)
    a.physical == b.physical && a.simulation == b.simulation &&
        a.solver == b.solver && a.hyper == b.hyper &&
        a.UDE == b.UDE && a.inversion == b.inversion
end

# Display setup
Base.show(io::IO, ::MIME"text/plain", params::Parameters) = Base.show(io, params)
function Base.show(io::IO, params::Parameters)
    label(s) = printstyled(io, rpad(s, 14); color = :light_black)
    sep() = printstyled(io, " · "; color = :light_black)
    field(s) = printstyled(io, s; color = :light_black)
    val(s) = print(io, s)
    hint(s) = printstyled(io, s; color = :light_black)
    check(b) = b ? "\e[32m✓\e[0m" : "\e[31m✗\e[0m"

    println(io, "Parameters")

    # Physical
    label("  physical")
    if isnothing(params.physical)
        hint("(nothing)")
    else
        p = params.physical
        field("ρ");
        print(io, " = ");
        val("$(p.ρ)")
        sep()
        field("A");
        print(io, " ∈ [");
        val("$(p.minA)");
        print(io, ", ");
        val("$(p.maxA)");
        print(io, "]")
        sep()
        field("C");
        print(io, " ∈ [");
        val("$(p.minC)");
        print(io, ", ");
        val("$(p.maxC)");
        print(io, "]")
    end
    println(io)

    # Simulation
    label("  simulation")
    if isnothing(params.simulation)
        hint("(nothing)")
    else
        s = params.simulation
        field("tspan");
        print(io, " = ");
        val("$(s.tspan)")
        sep()
        print(io, check(s.use_iceflow));
        field("iceflow")
        print(io, " ")
        print(io, check(s.use_MB));
        field("MB")
        print(io, " ")
        print(io, check(s.use_velocities));
        field("velocities")
    end
    println(io)

    # Solver
    label("  solver")
    if isnothing(params.solver)
        hint("(nothing)")
    else
        sv = params.solver
        val("$(nameof(typeof(sv.solver)))")
        sep()
        field("reltol");
        print(io, " = ");
        val("$(sv.reltol)")
        sep()
        field("maxiters");
        print(io, " = ");
        val("$(sv.maxiters)")
    end
    println(io)

    # Hyper
    label("  hyper")
    if isnothing(params.hyper)
        hint("(nothing)")
    else
        h = params.hyper
        field("epochs");
        print(io, " = ");
        val("$(h.epochs)")
        sep()
        field("batch_size");
        print(io, " = ");
        val("$(h.batch_size)")
        sep()
        if h.optimizer isa Vector
            opt_names = join([nameof(typeof(o)) for o in h.optimizer], ", ")
            val("[$(opt_names)]")
        else
            val("$(nameof(typeof(h.optimizer)))")
        end
    end
    println(io)

    # UDE
    label("  UDE")
    if isnothing(params.UDE)
        hint("(nothing)")
    else
        u = params.UDE
        field("target");
        print(io, " = ")
        isnothing(u.target) ? hint("(nothing)") : val(":$(u.target)")
        sep()
        val("\"$(u.optimization_method)\"")
        sep()
        isnothing(u.grad) ? hint("(nothing)") : val("$(nameof(typeof(u.grad)))")
        sep()
        field("loss");
        print(io, " = ");
        val("$(nameof(typeof(u.empirical_loss_function)))")
    end
    println(io)
end
