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

function generate_inputs(inputs::NamedTuple, simulation, glacier_idx, t)
    map(inputs) do input
        get_input(input, simulation, glacier_idx, t)
    end
end

"""
    WithInputs

Given tuple of `AbstractInput`s and a function `f`, return a callable struct 
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

function (f::WithInputs)(simulation, glacier_idx, t, θ)
    inputs = generate_inputs(f.inputs, simulation, glacier_idx, t)
    return f.f(inputs, θ)
end

function (f::WithInputs)(cache, simulation, glacier_idx, t, θ)
    inputs = generate_inputs(f.inputs, simulation, glacier_idx, t)
    return f.f(cache, inputs, θ)
end

_normalize_law_f(::Nothing, f) = f
_normalize_law_f(inputs, f) = WithInputs(inputs, f)

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
apply_law!(::AbstractLaw, cache, simulation, glacier_idx, t, θ) = throw(error("Concrete subtypes of AbstractLaw must implement apply_law!"))
init_cache(::AbstractLaw, simulation, glacier_idx, θ) = throw(error("Concrete subtypes of AbstractLaw must implement init_cache"))

"""
    Law{T}(;f!, init_cache, inputs = nothing, callback_freq = nothing)    

Defines a law applied to a model that mutates a cached internal state at each time step.

This struct encapsulates:
- a state update function `f!`,
- a cache initializer,
- and optionally a callback frequency (i.e. apply only every `callback_freq` time units).

Arguments:
- `f!::Function`: A function with the signature `f!(cache::T, simulation, glacier_idx, t, θ)` that updates the internal state (`cache`), return value is ignored.
  If `inputs` are provided, this function instead receives `f!(cache::T, inputs, θ)`
- `init_cache::Function`: Function with signature `init_cache(simulation, glacier_idx, θ)::T` to initialize the internal state.
- `callback_freq::Union{Nothing, AbstractFloat}`: Optional. If provided, the law becomes a callback law, applied every `callback_freq` time units.
- `inputs::Union{Nothing, Tuple{<:AbstractInput}}`: Optional. Allows automatic generation of inputs passed to `f!` at runtime.
"""
struct Law{CACHE_TYPE, F, INIT, FREQ} <: AbstractLaw
    f::F
    init_cache::INIT
    callback_freq::FREQ

    function Law{CACHE_TYPE}(f, init_cache, callback_freq) where {CACHE_TYPE}
        new{
            CACHE_TYPE,
            typeof(f),
            typeof(init_cache),
            typeof(callback_freq),
        }(
            f,
            init_cache,
            callback_freq,
        )
    end
end

ConstantLaw(value::T) where T = Law{T}(Returns(value), (_, _, _) -> (), nothing)

Law{T}(inputs, f, init_cache, callback_freq) where {T} = Law{T}(WithInputs(inputs, f), init_cache, callback_freq)
Law{T}(::Nothing, f, init_cache, callback_freq) where{T} = Law{T}(f, init_cache, callback_freq)
Law{T}(;f!, inputs = nothing, callback_freq = nothing, init_cache) where{T} = Law{T}(inputs, f!, init_cache, callback_freq)

apply_law!(law::Law, cache, simulation, glacier_idx, t, θ) = law.f(cache, simulation, glacier_idx, t, θ)
init_cache(law::Law, simulation, glacier_idx, θ) = law.init_cache(simulation, glacier_idx, θ)
cache_type(law::Law{CACHE_TYPE}) where {CACHE_TYPE} = CACHE_TYPE

is_callback_law(::Law{<:Any, <:Any, <:Any, Nothing}) = false
is_callback_law(::Law{<:Any, <:Any, <:Any, <:AbstractFloat}) = true

callback_freq(::Law{<:Any, <:Any, <:Any, Nothing}) = throw("This law dont have callback")
callback_freq(law::Law{<:Any, <:Any, <:Any, <:AbstractFloat}) = law.callback_freq

export Law
