export AbstractInput, get_input, generate_inputs, default_name
export AbstractPrepVJP

"""
    AbstractInput

Abstract type representing an input source for a `Law`.
Concrete subtypes must implement:

- `default_name(::ConcreteInput)`: returns the default field name (as a `Symbol`) under which this input will be accessed in the law's input `NamedTuple`.
- `get_input(::ConcreteInput, simulation, glacier_idx, t)`: retrieves the actual input value given the simulation context, glacier index, and time `t`.
"""
abstract type AbstractInput end

default_name(inp::AbstractInput) = throw(error("Concrete subtypes of AbstractInput must implement default_name. Please provide an implementation for $(typeof(inp))."))
get_input(inp::AbstractInput, simulation, glacier_idx, t) = throw(error("Concrete subtypes of AbstractInput must implement get_input. Please provide an implementation for $(typeof(inp))."))

function generate_inputs(inputs::NamedTuple, simulation, glacier_idx, t)
    map(inputs) do input
        get_input(input, simulation, glacier_idx, t)
    end
end

"""
    GenInputsAndApply{IN, F}

Given a tuple of `AbstractInput`s and a function `f`, returns a callable struct.
This struct `with_input_f` can be evaluated as a function that generates the inputs and then applies the function's law `with_input_f.f`.
It is for internal use only and it isn't exposed to the user.
"""
struct GenInputsAndApply{IN, F}
    inputs::IN
    f::F

    function GenInputsAndApply(inputs, f)
        inputs = _normalize_law_inputs(inputs)

        return new{typeof(inputs), typeof(f)}(inputs, f)
    end
end

_normalize_law_inputs(inputs::NamedTuple) = inputs
_normalize_law_inputs(inputs) = NamedTuple{default_name.(inputs)}(inputs)

function (with_input_f::GenInputsAndApply)(simulation, glacier_idx::Integer, t::Real, θ)
    inputs = generate_inputs(with_input_f.inputs, simulation, glacier_idx, t)
    return with_input_f.f(inputs, θ)
end

function (with_input_f::GenInputsAndApply)(cache, simulation, glacier_idx::Integer, t::Real, θ)
    inputs = generate_inputs(with_input_f.inputs, simulation, glacier_idx, t)
    return with_input_f.f(cache, inputs, θ)
end

"""
    AbstractPrepVJP

Abstract type representing the preparation of Vector-Jacobian Product (VJP) computations for the laws.
Subtypes of `AbstractPrepVJP` are used to handle any precomputations or setup required before evaluating VJPs, such as configuring automatic differentiation backends or precompiling code.
This type provides a flexible interface for implementing custom VJP preparation strategies for different laws.
It is for internal use only and not exposed to the user.
"""
abstract type AbstractPrepVJP end

function (with_input_f::GenInputsAndApply)(vjpPrep::Union{AbstractPrepVJP, Nothing}, simulation, glacier_idx::Integer, t::Real, θ)
    inputs = generate_inputs(with_input_f.inputs, simulation, glacier_idx, t)
    return with_input_f.f(vjpPrep, inputs, θ)
end

function (with_input_f::GenInputsAndApply)(cache, vjpPrep::Union{AbstractPrepVJP, Nothing}, simulation, glacier_idx::Integer, t::Real, θ)
    inputs = generate_inputs(with_input_f.inputs, simulation, glacier_idx, t)
    return with_input_f.f(cache, vjpPrep, inputs, θ)
end
