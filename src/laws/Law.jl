export NullLaw, AbstractLaw, Law, ConstantLaw
export init_cache, cache_type, is_callback_law, callback_freq, build_affect, apply_law!, apply_law_in_model, apply_all_non_callback_laws!, apply_all_callback_laws!, precompute_all_VJPs_laws!
export AbstractInput, get_input, inputs, generate_inputs, default_name
export law_VJP_input, law_VJP_θ, prepare_vjp_law, precompute_law_VJP, is_precomputable_law_VJP
export CustomVJP, DIVJP
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

function (with_input_f::GenInputsAndApply)(vjpPrep::AbstractPrepVJP, simulation, glacier_idx::Integer, t::Real, θ)
    inputs = generate_inputs(with_input_f.inputs, simulation, glacier_idx, t)
    return with_input_f.f(vjpPrep, inputs, θ)
end

function (with_input_f::GenInputsAndApply)(cache, vjpPrep::AbstractPrepVJP, simulation, glacier_idx::Integer, t::Real, θ)
    inputs = generate_inputs(with_input_f.inputs, simulation, glacier_idx, t)
    return with_input_f.f(cache, vjpPrep, inputs, θ)
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
    emptyVJP(cache, simulation, glacier_idx, t, θ)
    emptyVJPWithInputs(cache, inputs, θ)

Function that defines an empty `!`-style function for the VJPs.
The two methods define the two possible signatures for the function that updates the cache in-place.
It is for internal use only and it isn't exposed to the user.
"""
emptyVJP(cache, simulation, glacier_idx, t, θ) = nothing
emptyVJPWithInputs(cache, inputs, θ) = nothing

"""
    emptyPrepVJP(cache, vjpPrep, simulation, glacier_idx, t, θ)
    emptyPrepVJPWithInputs(cache, vjpPrep, inputs, θ)

Function that defines an empty `!`-style function for the preparation of the VJPs.
The two methods define the two possible signatures for the function that updates the cache in-place.
Trying to apply this function will yield an error.
It is for internal use only and it isn't exposed to the user.
"""
emptyPrepVJP(cache, vjpPrep, simulation, glacier_idx, t, θ) = throw("This VJP preparation has not been defined.")
emptyPrepVJPWithInputs(cache, vjpPrep, inputs, θ) = throw("This VJP preparation has not been defined.")

"""
    ∂law∂θ!()

This function serves as a placeholder and should be replaced by other implementations in ODINN.
This implementation throws an error.
It is for internal use only and is not exposed to the user.
"""
∂law∂θ!() = throw("Function ∂law∂θ! not implemented. Please provide an implementation.")

"""
    ∂law∂inp!()

This function serves as a placeholder and should be replaced by other implementations in ODINN.
This implementation throws an error.
It is for internal use only and is not exposed to the user.
"""
∂law∂inp!() = throw("Function ∂law∂inp! not implemented. Please provide an implementation.")

"""
    prepare_vjp_law(simulation, law::AbstractLaw, law_cache, θ, glacier_idx)

Function used to prepare the VJPs at the initialization of the model cache.
It is used for example to compile VJPs of the laws to be differentiated using DifferentiationInterface.jl.
"""
prepare_vjp_law(simulation, law::AbstractLaw, law_cache, θ, glacier_idx) = throw("This prepare function has not been defined.")

"""
    VJPType

Abstract type representing the mode of Vector-Jacobian Product (VJP) computation used in a law.
Subtypes of `VJPType` define how the VJP is evaluated (e.g., custom or using DifferentiationInterface.jl).
"""
abstract type VJPType end

"""
    CustomVJP <: VJPType

Indicates that a law uses a custom-defined function for VJP computation.
This is used when the VJP is provided manually rather than computed automatically.
"""
struct CustomVJP <: VJPType end

"""
    DIVJP <: VJPType

Indicates that a law uses the default VJP computation provided by DifferentiationInterface.jl.
This is used when no custom VJP function is provided, and the Vector-Jacobian Product (VJP) is computed automatically using DifferentiationInterface.jl.
"""
struct DIVJP <: VJPType end

"""
    Law{T}(; f!, init_cache, callback_freq=nothing, inputs=nothing)

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
- `callback_freq::Union{Nothing, Real}`: Optional. If provided, the law is treated as a callback law and is only applied every `callback_freq` time units. If `callback_freq` is set to zero, then the law is applied only once at the beginning of the simulation.
- `inputs::Union{Nothing, Tuple{<:AbstractInput}}`: Optional. Provides automatically generated inputs passed to `f!` at runtime.

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
struct Law{
    CACHE_TYPE,
    F, FVJPINP, FVJPθ,
    INIT,
    FREQ,
    PFVJP,
    VJPType,
} <: AbstractLaw
    f::F
    f_VJP_input::FVJPINP
    f_VJP_θ::FVJPθ
    init_cache::INIT
    callback_freq::FREQ
    p_VJP::PFVJP
    VJPType::VJPType

    function Law{CACHE_TYPE}(
        f, f_VJP_input, f_VJP_θ, init_cache, callback_freq, p_VJP, vjpType
    ) where {CACHE_TYPE}
        new{
            CACHE_TYPE,
            typeof(f),
            typeof(f_VJP_input),
            typeof(f_VJP_θ),
            typeof(init_cache),
            typeof(callback_freq),
            typeof(p_VJP),
            typeof(vjpType)
        }(
            f,
            f_VJP_input,
            f_VJP_θ,
            init_cache,
            callback_freq,
            p_VJP,
            vjpType,
        )
    end
end


### Declaration of the law with semantic inputs (cf AbstractInput) ###
Law{T}( # No custom VJP
    inputs::Union{<: NamedTuple, <: Tuple}, f::Function, init_cache, callback_freq
) where {T} = Law{T}(
    GenInputsAndApply(inputs, f),
    GenInputsAndApply(inputs, emptyVJPWithInputs),
    GenInputsAndApply(inputs, emptyVJPWithInputs),
    init_cache, callback_freq,
    GenInputsAndApply(inputs, emptyPrepVJPWithInputs),
    DIVJP(),
)
Law{T}( # With VJP computed on-the-fly
    inputs::Union{<: NamedTuple, <: Tuple},
    f::Function, f_VJP_input::Function, f_VJP_θ::Function,
    init_cache, callback_freq,
) where {T} = Law{T}(
    GenInputsAndApply(inputs, f),
    GenInputsAndApply(inputs, f_VJP_input),
    GenInputsAndApply(inputs, f_VJP_θ),
    init_cache, callback_freq,
    GenInputsAndApply(inputs, emptyPrepVJPWithInputs),
    CustomVJP(),
)
Law{T}( # With precomputed VJP and VJP interpolated on-the-fly
    inputs::Union{<: NamedTuple, <: Tuple},
    f::Function, f_VJP_input::Function, f_VJP_θ::Function,
    init_cache, callback_freq,
    p_VJP::Union{Function, DIVJP},
) where {T} = Law{T}(
    GenInputsAndApply(inputs, f),
    GenInputsAndApply(inputs, f_VJP_input),
    GenInputsAndApply(inputs, f_VJP_θ),
    init_cache, callback_freq,
    GenInputsAndApply(inputs, p_VJP),
    CustomVJP(),
)
Law{T}( # With precomputed VJP and no VJP interpolation on-the-fly
    inputs::Union{<: NamedTuple, <: Tuple},
    f::Function, f_VJP_input::Nothing, f_VJP_θ::Nothing,
    init_cache, callback_freq,
    p_VJP::Union{Function, DIVJP},
) where {T} = Law{T}(
    GenInputsAndApply(inputs, f),
    GenInputsAndApply(inputs, emptyVJPWithInputs),
    GenInputsAndApply(inputs, emptyVJPWithInputs),
    init_cache, callback_freq,
    GenInputsAndApply(inputs, p_VJP),
    CustomVJP(),
)
Law{T}( # Binding for the case with no custom VJP
    inputs::Union{<: NamedTuple, <: Tuple},
    f::Function, f_VJP_input::Nothing, f_VJP_θ::Nothing,
    init_cache, callback_freq,
    p_VJP::Nothing,
) where {T} = Law{T}(inputs, f, init_cache, callback_freq)
Law{T}( # Binding for the case with VJP computed on-the-fly
    inputs::Union{<: NamedTuple, <: Tuple},
    f::Function, f_VJP_input::Function, f_VJP_θ::Function,
    init_cache, callback_freq,
    p_VJP::Nothing,
) where {T} = Law{T}(inputs, f, f_VJP_input, f_VJP_θ, init_cache, callback_freq)


### Declaration of the law with an affect that directly retrieves the inputs ###
Law{T}( # No custom VJP
    ::Nothing, f::Function, init_cache, callback_freq
) where{T} = Law{T}(
    f, emptyVJP, emptyVJP,
    init_cache, callback_freq,
    emptyPrepVJP,
    DIVJP(),
)
Law{T}( # With VJP computed on-the-fly
    ::Nothing,
    f::Function, f_VJP_input::Function, f_VJP_θ::Function,
    init_cache, callback_freq,
) where{T} = Law{T}(
    f, f_VJP_input, f_VJP_θ,
    init_cache, callback_freq,
    emptyPrepVJP,
    CustomVJP(),
)
Law{T}( # With precomputed VJP and VJP interpolated on-the-fly
    ::Nothing,
    f::Function, f_VJP_input::Function, f_VJP_θ::Function,
    init_cache, callback_freq,
    p_VJP::Function,
) where{T} = Law{T}(
    f, f_VJP_input, f_VJP_θ,
    init_cache, callback_freq,
    p_VJP,
    CustomVJP(),
)
Law{T}( # With precomputed VJP and no VJP interpolation on-the-fly
    ::Nothing,
    f::Function, f_VJP_input::Nothing, f_VJP_θ::Nothing,
    init_cache, callback_freq,
    p_VJP::Function,
) where{T} = Law{T}(
    f, emptyVJP, emptyVJP,
    init_cache, callback_freq,
    p_VJP,
    CustomVJP(),
)
Law{T}( # Binding for the case with no custom VJP
    ::Nothing,
    f::Function, f_VJP_input::Nothing, f_VJP_θ::Nothing,
    init_cache, callback_freq,
    p_VJP::Nothing,
) where {T} = Law{T}(nothing, f, init_cache, callback_freq)
Law{T}( # Binding for the case with VJP computed on-the-fly
    ::Nothing,
    f::Function, f_VJP_input::Function, f_VJP_θ::Function,
    init_cache, callback_freq,
    p_VJP::Nothing,
) where {T} = Law{T}(f, f_VJP_input, f_VJP_θ, init_cache, callback_freq)

### Wrapper that calls the declaration of the law with semantic inputs (useful to use default arguments) ###
Law{T}(; # This definition relies on the bindings above
    f!, f_VJP_input! = nothing, f_VJP_θ! = nothing,
    inputs = nothing, callback_freq = nothing,
    init_cache,
    p_VJP! = nothing,
) where{T} = Law{T}(
    inputs,
    f!, f_VJP_input!, f_VJP_θ!,
    init_cache, callback_freq,
    p_VJP!,
)

apply_law!(law::Law, cache, simulation, glacier_idx, t, θ) = law.f(cache, simulation, glacier_idx, t, θ)
law_VJP_input(law::Law, cache, simulation, glacier_idx, t, θ) = law.f_VJP_input(cache, simulation, glacier_idx, t, θ)
law_VJP_θ(law::Law, cache, simulation, glacier_idx, t, θ) = law.f_VJP_θ(cache, simulation, glacier_idx, t, θ)
precompute_law_VJP(law::Law, cache, vjpsPrepLaw::AbstractPrepVJP, simulation, glacier_idx, t, θ) = law.p_VJP(cache, vjpsPrepLaw, simulation, glacier_idx, t, θ)
function precompute_law_VJP(law::Law{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:GenInputsAndApply{<:Any, DIVJP}, <:Any}, cache, vjpsPrepLaw::AbstractPrepVJP, simulation, glacier_idx, t, θ)
    backend = simulation.parameters.UDE.grad.VJP_method.regressorADBackend
    inputs = generate_inputs(law.p_VJP.inputs, simulation, glacier_idx, t)
    ∂θ = ∂law∂θ!(backend, vjpsPrepLaw, inputs, θ)
    cache.vjp_θ .= ∂θ
    ∂inp = ∂law∂inp!(backend, vjpsPrepLaw, inputs, θ)
    cache.vjp_inp .= ∂inp
end

init_cache(law::Law, simulation, glacier_idx, θ) = law.init_cache(simulation, glacier_idx, θ)
cache_type(law::Law{CACHE_TYPE}) where {CACHE_TYPE} = CACHE_TYPE

is_callback_law(::Law{<:Any, <:Any, <:Any, <:Any, <:Any, Nothing, <:Any, <:Any}) = false
is_callback_law(::Law{<:Any, <:Any, <:Any, <:Any, <:Any, <:Real, <:Any, <:Any}) = true
is_precomputable_law_VJP(law::Law{CACHE_TYPE, <:Any, <:Any, <:Any, <:Any, <:Any, <:GenInputsAndApply, <:Any}) where {CACHE_TYPE} = !isa(law.p_VJP.f, typeof(emptyPrepVJPWithInputs))
is_precomputable_law_VJP(law::Law) = !(isa(law.p_VJP, typeof(emptyPrepVJP)) || isnothing(law.p_VJP))

callback_freq(::Law{<:Any, <:Any, <:Any, <:Any, <:Any, Nothing, <:Any}) = throw("This law does not have callback")
callback_freq(law::Law{<:Any, <:Any, <:Any, <:Any, <:Any, <:Real, <:Any}) = law.callback_freq

# Define whether inputs are provided or not through GenInputsAndApply
inputs(law::Law{CACHE_TYPE, <: GenInputsAndApply}) where {CACHE_TYPE} = law.f.inputs
inputs(law::Law) = throw("Inputs are not defined.")
apply_law_in_model(law::Law) = !is_callback_law(law)


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
law_VJP_input(law::ConstantLaw, cache, simulation, glacier_idx, t, θ) = nothing
law_VJP_θ(law::ConstantLaw, cache, simulation, glacier_idx, t, θ) = nothing
precompute_law_VJP(law::ConstantLaw, cache, vjpsPrepLaw::AbstractPrepVJP, simulation, glacier_idx, t, θ) = nothing
init_cache(law::ConstantLaw, simulation, glacier_idx, θ) = law.init_cache(simulation, glacier_idx, θ)
cache_type(law::ConstantLaw{CACHE_TYPE}) where {CACHE_TYPE} = CACHE_TYPE
is_callback_law(::ConstantLaw) = false
is_precomputable_law_VJP(law::ConstantLaw) = false
callback_freq(::ConstantLaw) = throw("ConstantLaw doesn't have callback")
inputs(law::ConstantLaw) = (;)
apply_law_in_model(law::ConstantLaw) = false


"""
    NullLaw <: AbstractLaw

This struct represents a law that is not used in the iceflow model.
"""
struct NullLaw <: AbstractLaw end

apply_law!(law::NullLaw, cache, simulation, glacier_idx, t, θ) = nothing
law_VJP_input(law::NullLaw, cache, simulation, glacier_idx, t, θ) = nothing
law_VJP_θ(law::NullLaw, cache, simulation, glacier_idx, t, θ) = nothing
precompute_law_VJP(law::NullLaw, cache, vjpsPrepLaw::AbstractPrepVJP, simulation, glacier_idx, t, θ) = nothing
init_cache(law::NullLaw, simulation, glacier_idx, θ) = FloatCacheNoVJP(zeros(Sleipnir.Float))
cache_type(law::NullLaw) = FloatCacheNoVJP
is_callback_law(::NullLaw) = false
is_precomputable_law_VJP(law::NullLaw) = false
callback_freq(::NullLaw) = throw("NullLaw doesn't have callback")
inputs(law::NullLaw) = (;)
apply_law_in_model(law::NullLaw) = false

"""
    apply_all_non_callback_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ)

This function is a placeholder and must be implemented for your custom model type.

It is intended to apply all non-callback laws in the simulation to update the model cache for a given glacier at time `t` with parameters `θ`. By default, calling this function will throw an error to indicate that the user should provide their own implementation tailored to their model.

# Arguments
- `model::AbstractModel`: The model instance.
- `cache`: The cache object storing state variables.
- `simulation`: The simulation context.
- `glacier_idx`: Index identifying the glacier.
- `t`: The current simulation time.
- `θ`: The parameter vector.

# Throws
- Always throws an error: `"This function should not be called. Implement apply_all_non_callback_laws! for your own model."`
"""
apply_all_non_callback_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ) = throw("This function should not be called. Implement apply_all_non_callback_laws! for your own model.")

"""
    precompute_all_VJPs_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ)

This function is a placeholder and must be implemented for your custom model type.

It is intended to precompute the VJPs for all the laws that are used in a model. By default, calling this function will throw an error to indicate that the user should provide their own implementation tailored to their model.

# Arguments
- `model::AbstractModel`: The model instance.
- `cache`: The cache object storing state variables.
- `simulation`: The simulation context.
- `glacier_idx`: Index identifying the glacier.
- `t`: The current simulation time.
- `θ`: The parameter vector.

# Throws
- Always throws an error: `"This function should not be called. Implement precompute_all_VJPs_laws! for your own model."`
"""
precompute_all_VJPs_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ) = throw("This function should not be called. Implement precompute_all_VJPs_laws! for your own model.")

"""
    apply_all_callback_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ)

This function is a placeholder and must be implemented for your custom model type.

It is intended to apply all callback laws in the simulation to update the model cache for a given glacier at time `t` with parameters `θ`. By default, calling this function will throw an error to indicate that the user should provide their own implementation tailored to their model.

# Arguments
- `model::AbstractModel`: The model instance.
- `cache`: The cache object storing state variables.
- `simulation`: The simulation context.
- `glacier_idx`: Index identifying the glacier.
- `t`: The current simulation time.
- `θ`: The parameter vector.

# Throws
- Always throws an error: `"This function should not be called. Implement apply_all_callback_laws! for your own model."`
"""
apply_all_callback_laws!(model::AbstractModel, cache, simulation, glacier_idx, t, θ) = throw("This function should not be called. Implement apply_all_callback_laws! for your own model.")


# Display setup
repr(law::Law{CACHE_TYPE, <: GenInputsAndApply}) where {CACHE_TYPE} = "$(keys(inputs(law))) -> $(cache_type(law))"
repr(law::Law) = "(simulation, t) -> $(cache_type(law))"
repr(law::ConstantLaw) = "ConstantLaw -> $(cache_type(law))"
repr(law::NullLaw) = "NullLaw"
Base.show(io::IO, type::MIME"text/plain", law::AbstractLaw) = Base.show(io, law)
function Base.show(io::IO, law::AbstractLaw)
    print(io, repr(law))
end
