"""
    AbstractInput

Abstract type representing an input source for a `Law`.  
Concrete subtypes must implement:

- `default_name(::ConcreteInput)`: returns the default field name (as a `Symbol`) under which this input will be accessed in the law's input `NamedTuple`.
- `get_input(::ConcreteInput, simulation, glacier_idx, t)`: retrieves the actual input value given the simulation context, glacier index, and time `t`.
"""
abstract type AbstractInput end

default_name(::AbstractInput) = throw(error("Concrete subtypes of AbstractInput must implement default_name"))
get_input(::AbstractInput, simulation, glacier_idx, t) = throw(error("Concrete subtypes of AbstractInput must implement get_input"))

"""
    AbstractLaw

Abstract type representing a synthetic law.
Currently it's only used for testing by making easier to create dumb laws, but in the future it may be cleaner to use different concrete type of laws
(for example `CallbackLaw`, `ContinuousLaw`, or `LearnableLaw`)

Concrete subtypes must implement:
- `apply_law(::ConcreteLaw, simulation, glacier_idx, t, θ)`
"""
abstract type AbstractLaw end

"""
    Law(; inputs::Tuple|NamedTuple, f::Function, callback_freq)

A `Law` can represent either a synthetic function or a learnable model (e.g., neural network) that computes a parameter based on selected input variables.

# Arguments
- `inputs`: a `Tuple` or `NamedTuple` of `AbstractInput` subtypes describing the data dependencies of the law or learner.
- `f`: the function defining the law or learner logic. It should accept two arguments:
    1. A `NamedTuple` of input values (with the same field names as the `inputs`).
    2. A parameter vector `θ` (e.g., model parameters).
- `callback_freq`: a `Float64` specifying how often (in simulation time units) the law should be computed. Use `nothing` or `0.0` for continuous application.

# Example

```julia
law = Law(
    inputs = (Input1(), Input2()),
    f = (inputs, θ) -> inputs.input1 + θ[1] * inputs.input2,
    callback = 1.0
)
```
"""
struct Law{IN, F, FREQ <: Union{Nothing, Float64}} <: AbstractLaw
    inputs::IN
    f::F
    callback_freq::FREQ

    function Law(; inputs, f, callback_freq = nothing)
        # internally we always store inputs as NamedTuple
        inputs = _normalize_law_inputs(inputs)

        new{typeof(inputs), typeof(f), typeof(callback_freq)}(inputs, f, callback_freq)
    end
end

_normalize_law_inputs(inputs::NamedTuple) = inputs
_normalize_law_inputs(inputs) = NamedTuple{default_name.(inputs)}(inputs)

function generate_inputs(inputs::NamedTuple, simulation, glacier_idx, t)
    map(inputs) do input
        get_input(input, simulation, glacier_idx, t)
    end
end

function apply_law(law::Law, simulation, glacier_idx, t, θ)
    inputs = generate_inputs(law.inputs, simulation, glacier_idx, t)
    law.f(inputs, θ)
end
