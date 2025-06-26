export NullLaw, AbstractLaw, Law, ConstantLaw, init_cache, cache_type, is_differentiable, is_callback_law, callback_freq, build_affect, apply_law!, get_input, AbstractInput, generate_inputs

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
    WithInputs{IN, F}

Given a tuple of `AbstractInput`s and a function `f`, returns a callable struct.
This struct `with_input_f` can be evaluated as a function that generates the inputs and then applies the function's law `with_input_f.f`.
It is for internal use only and is not exposed to the user.
"""
struct WithInputs{IN, F}
    inputs::IN
    f::F

    function WithInputs(inputs, f)
        inputs = _normalize_law_inputs(inputs)

        return new{typeof(inputs), typeof(f)}(inputs, f)
    end
end

_normalize_law_inputs(inputs::NamedTuple) = inputs
_normalize_law_inputs(inputs) = NamedTuple{default_name.(inputs)}(inputs)

function (with_input_f::WithInputs)(simulation, glacier_idx, t, θ)
    inputs = generate_inputs(with_input_f.inputs, simulation, glacier_idx, t)
    return with_input_f.f(inputs, θ)
end

function (with_input_f::WithInputs)(cache, simulation, glacier_idx, t, θ)
    inputs = generate_inputs(with_input_f.inputs, simulation, glacier_idx, t)
    return with_input_f.f(cache, inputs, θ)
end

"""
    AbstractLaw

Abstract type representing a synthetic law.
Currently it's only used for testing by making easier to create dumb laws, but in the future it may be cleaner to use different concrete type of laws
(for example `CallbackLaw`, `ContinuousLaw`, or `LearnableLaw`)

Concrete subtypes must implement:
- `apply_law!(::ConcreteLaw, state, simulation, glacier_idx, t, θ)`
- `init_cache(::ConcreteLaw, glacier, glacier_idx)::`
"""
abstract type AbstractLaw end
apply_law!(inp::AbstractLaw, cache, simulation, glacier_idx, t, θ) = throw(error("Concrete subtypes of AbstractLaw must implement apply_law!. Please provide an implementation for $(typeof(inp))."))
init_cache(inp::AbstractLaw, simulation, glacier_idx, θ) = throw(error("Concrete subtypes of AbstractLaw must implement init_cache. Please provide an implementation for $(typeof(inp))."))

"""
    build_affect(law::AbstractLaw, cache, glacier_idx, θ)

Return a `!`-style function suitable for use in a callback, which applies the given `law`
to update the `cache` for a specific glacier and parameters `θ`, using the simulation time.
"""
function build_affect(law::AbstractLaw, cache, glacier_idx, θ)
    # The let block make sure that every variable are type stable
    return let law = law, cache = cache, glacier_idx = glacier_idx, θ = θ
        function affect!(integrator)
            simulation = integrator.p
            t = integrator.t

            apply_law!(law, cache, simulation, glacier_idx, t, θ)
        end
    end
end

"""
    Law{T}(; f!, init_cache, callback_freq=nothing, inputs=nothing, is_differentiable=false)

Defines a physical or empirical law applied to a glacier model that mutates an internal state `T` at each simulation time step.

!!! warning
    The type `T` must be *mutable*, since `f!` is expected to update `cache::T` in-place.
    Using an immutable type (like `Float64`) will silently fail or raise an error.

```
# ❌ Will not work: Float64 is immutable, so cache .= ... has no effect
Law{Float64}(;
    f! = (cache, _, _, t, θ) -> cache = θ.scale * sin(2π * t + θ.shift),
    init_cache = (_, _, _) -> 0.0,
)

# ✅ Correct: using a 0-dimensional array allows in-place mutation
Law{Array{Float64, 0}}(;
    f! = (cache, _, _, t, θ) -> cache .= θ.scale * sin(2π * t + θ.shift) + θ.bias,
    init_cache = (_, _, _) -> zeros(),
)
```

# Arguments

- `f!::Function`: A function with signature `f!(cache::T, simulation, glacier_idx, t, θ)` that updates the internal state.
  If `inputs` are provided, the function instead takes the form `f!(cache::T, inputs, θ)`.
- `init_cache::Function`: A function `init_cache(simulation, glacier_idx, θ)::T` that initializes the internal state for a given glacier.
- `callback_freq::Union{Nothing, AbstractFloat}`: Optional. If provided, the law is treated as a callback law and is only applied every `callback_freq` time units.
- `inputs::Union{Nothing, Tuple{<:AbstractInput}}`: Optional. Provides automatically generated inputs passed to `f!` at runtime.
- `is_differentiable::Bool`: Optional. Whether the law can be differentiated within ODINN or not.
  For a law to be differentiable within ODINN, its inputs must be carefully handled.
  For the moment only a subset of laws are supported and the users should not use this
  parameter unless they know what they are doing. Trying to use unsupported laws for
  inversions in ODINN is highly discouraged as it may result in incorrect gradients.
  In the future we plan to fully support arbitrary laws once ODINN can be fully
  differentiated using Enzyme.

# Type Parameters

- `T`: The type of the internal state. Must be specified manually and should match the return type of `init_cache`.

# Examples

```
# A law applied at every timestep, storing a scalar value
Law{Array{Float64, 0}}(;
    f! = (cache, _, _, t, θ) -> cache .= θ.scale * sin(2π * t + θ.shift) + θ.bias,
    init_cache = (_, _, _) -> zeros(),
)

# A callback law applied once per month (assuming time in years)
Law{Array{Float64, 0}}(;
    f! = (cache, _, _, t, θ) -> cache .= θ.scale * sin(2π * t + θ.shift) + θ.bias,
    init_cache = (_, _, _) -> zeros(),
    callback_freq = 1 / 12,
)

# TODO: create a meaningful example of Law with glacier-shaped matrix
```
"""
struct Law{CACHE_TYPE, F, INIT, FREQ} <: AbstractLaw
    f::F
    init_cache::INIT
    callback_freq::FREQ
    is_differentiable::Bool

    function Law{CACHE_TYPE}(f, init_cache, callback_freq; is_differentiable::Bool = false) where {CACHE_TYPE}
        new{
            CACHE_TYPE,
            typeof(f),
            typeof(init_cache),
            typeof(callback_freq),
        }(
            f,
            init_cache,
            callback_freq,
            is_differentiable,
        )
    end
end

Law{T}(inputs, f, init_cache, callback_freq; is_differentiable::Bool = false) where {T} = Law{T}(WithInputs(inputs, f), init_cache, callback_freq; is_differentiable=is_differentiable)
Law{T}(::Nothing, f, init_cache, callback_freq; is_differentiable::Bool = false) where{T} = Law{T}(f, init_cache, callback_freq; is_differentiable=is_differentiable)
Law{T}(;f!, inputs = nothing, callback_freq = nothing, init_cache, is_differentiable::Bool = false) where{T} = Law{T}(inputs, f!, init_cache, callback_freq; is_differentiable=is_differentiable)

apply_law!(law::Law, cache, simulation, glacier_idx, t, θ) = law.f(cache, simulation, glacier_idx, t, θ)
init_cache(law::Law, simulation, glacier_idx, θ) = law.init_cache(simulation, glacier_idx, θ)
cache_type(law::Law{CACHE_TYPE}) where {CACHE_TYPE} = CACHE_TYPE
is_differentiable(law::Law) = law.is_differentiable

is_callback_law(::Law{<:Any, <:Any, <:Any, Nothing}) = false
is_callback_law(::Law{<:Any, <:Any, <:Any, <:AbstractFloat}) = true

callback_freq(::Law{<:Any, <:Any, <:Any, Nothing}) = throw("This law does not have callback")
callback_freq(law::Law{<:Any, <:Any, <:Any, <:AbstractFloat}) = law.callback_freq


"""
    ConstantLaw{T}(init_cache)

Creates a constant law of type `Law{T}` that holds a fixed value for the entire simulation.

This is useful to inject glacier-specific or global constants into the simulation without modifying them over time.
The update function is a no-op, and only the `init_cache` function matters.

# Arguments

- `init_cache::Function`: A function `init_cache(simulation, glacier_idx, θ)::T` that provides the constant value.

# Type Parameters

- `T`: The type of the internal state. Must be specified manually and should match the return type of `init_cache`.

# Examples

```
# Same value for all glaciers
n_law = ConstantLaw{Float64}(Returns(4.))

# Value depending on the glacier
n_law = ConstantLaw{Float64}((sim, i, θ) -> sim.glaciers[i].n)

# Learned value
n_law = ConstantLaw{Float64}((sim, i, θ) -> θ.n)
```
"""
struct ConstantLaw{CACHE_TYPE, INIT} <: AbstractLaw
    init_cache::INIT

    function ConstantLaw{CACHE_TYPE}(init_cache) where CACHE_TYPE
        return new{CACHE_TYPE, typeof(init_cache)}(init_cache)
    end
end

apply_law!(law::ConstantLaw, cache, simulation, glacier_idx, t, θ) = nothing
init_cache(law::ConstantLaw, simulation, glacier_idx, θ) = law.init_cache(simulation, glacier_idx, θ)
cache_type(law::ConstantLaw{CACHE_TYPE}) where {CACHE_TYPE} = CACHE_TYPE
is_differentiable(law::ConstantLaw) = true
is_callback_law(::ConstantLaw) = false
callback_freq(::ConstantLaw) = throw("ConstantLaw doesn't have callback")


"""
    NullLaw <: AbstractLaw

This struct represents a law that is not used in the iceflow model.
"""
struct NullLaw <: AbstractLaw end

apply_law!(law::NullLaw, cache, simulation, glacier_idx, t, θ) = nothing
init_cache(law::NullLaw, simulation, glacier_idx, θ) = zeros(Sleipnir.Float)
cache_type(law::NullLaw) = Array{Sleipnir.Float, 0}
is_differentiable(law::NullLaw) = false
is_callback_law(::NullLaw) = false
callback_freq(::NullLaw) = throw("NullLaw doesn't have callback")
