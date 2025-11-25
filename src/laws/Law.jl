export NullLaw, Law, ConstantLaw

const CB_FREQ_TYPE = Union{Nothing, Real, Period, Month}

convert_callbackfreq(callback_freq::Union{Nothing, Real}) = callback_freq
convert_callbackfreq(callback_freq::Union{Period, Month}) = yearfrac(callback_freq)

"""
    Law{T}(;
        inputs = nothing,
        f!, f_VJP_input! = nothing, f_VJP_θ! = nothing,
        init_cache, callback_freq = nothing,
        p_VJP! = nothing,
        max_value = NaN, min_value = NaN, name = :unknown,
    ) where{T}

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
  - `callback_freq::Union{Nothing, Real, Period, Month}`: Optional. If provided, the law is treated as a callback law and is only applied every `callback_freq` time units.
    If `callback_freq` is set to zero, then the law is applied only once at the beginning of the simulation.
    If `callback_freq` is set to `nothing` (default), then the law is applied at every iteration.
    If `callback_freq` is provided as a `Period` or `Month`, it is converted to a float value in a yearly basis.
  - `f_VJP_input!`: A function with signature `(cache::T, simulation, glacier_idx, t, θ)` that updates `cache.vjp_inp` which is the VJP with respect to the inputs.
  - `f_VJP_θ!`: A function with signature `(cache::T, simulation, glacier_idx, t, θ)` that updates `cache.vjp_θ` which is the VJP with respect to the parameters θ.
  - `p_VJP!`: A function with signature `(cache::T, vjpsPrepLaw, simulation, glacier_idx, t, θ)` that performs the precomputation of the VJPs.
  - `inputs::Union{Nothing, Tuple{<:AbstractInput}}`: Optional. Provides automatically generated inputs passed to `f!` at runtime.
  - `max_value::Float64`: Optional. The maximum value that the law can take, used for plotting and capping the function output.
  - `min_value::Float64`: Optional. The minimum value that the law can take, used for plotting and capping the function output.
  - `name::Symbol`: A name for the law, used for identification and plotting.

# Type Parameters

  - `T`: The type of the internal state. Must be specified manually and should match the return type of `init_cache`.

# Notes

  - Refer to the tutorials in the documentation for a complete description of the VJP options.

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
```
"""
struct Law{
    CACHE_TYPE,
    F, FVJPINP, FVJPθ,
    INIT,
    FREQ,
    PFVJP,
    VJPType
} <: AbstractLaw
    f::F
    f_VJP_input::FVJPINP
    f_VJP_θ::FVJPθ
    init_cache::INIT
    callback_freq::FREQ
    p_VJP::PFVJP
    VJPType::VJPType
    max_value::Float64
    min_value::Float64
    name::Symbol

    # Inner constructor that will be called by all the other dispatch
    function Law{CACHE_TYPE}(
            f, f_VJP_input, f_VJP_θ, init_cache, callback_freq, p_VJP, vjpType,
            max_value, min_value, name
    ) where {CACHE_TYPE}
        converted_callback_freq = convert_callbackfreq(callback_freq)
        new{
            CACHE_TYPE,
            typeof(f),
            typeof(f_VJP_input),
            typeof(f_VJP_θ),
            typeof(init_cache),
            typeof(converted_callback_freq),
            typeof(p_VJP),
            typeof(vjpType)
        }(
            f,
            f_VJP_input,
            f_VJP_θ,
            init_cache,
            converted_callback_freq,
            p_VJP,
            vjpType,
            max_value,
            min_value,
            name
        )
    end
end

### A. Declaration of the law with semantic inputs (cf AbstractInput) ###
function Law{T}( # 1. No custom VJP
        inputs::Union{<: NamedTuple, <: Tuple},
        f::Function, f_VJP_input::Nothing, f_VJP_θ::Nothing,
        init_cache, callback_freq::CB_FREQ_TYPE,
        p_VJP::Nothing,
        max_value::Real, min_value::Real, name::Symbol
) where {T}
    Law{T}(
        GenInputsAndApply(inputs, f),
        GenInputsAndApply(inputs, emptyVJPWithInputs),
        GenInputsAndApply(inputs, emptyVJPWithInputs),
        init_cache, callback_freq,
        GenInputsAndApply(inputs, emptyPrepVJPWithInputs),
        DIVJP(),
        max_value, min_value, name
    )
end
function Law{T}( # 2. With VJP computed on-the-fly
        inputs::Union{<: NamedTuple, <: Tuple},
        f::Function, f_VJP_input::Function, f_VJP_θ::Function,
        init_cache, callback_freq::CB_FREQ_TYPE,
        p_VJP::Nothing,
        max_value::Real, min_value::Real, name::Symbol
) where {T}
    Law{T}(
        GenInputsAndApply(inputs, f),
        GenInputsAndApply(inputs, f_VJP_input),
        GenInputsAndApply(inputs, f_VJP_θ),
        init_cache, callback_freq,
        GenInputsAndApply(inputs, emptyPrepVJPWithInputs),
        CustomVJP(),
        max_value, min_value, name
    )
end
function Law{T}( # 3. With precomputed VJP and possibly some VJP computations on-the-fly
        inputs::Union{<: NamedTuple, <: Tuple},
        f::Function, f_VJP_input::Function, f_VJP_θ::Function,
        init_cache, callback_freq::CB_FREQ_TYPE,
        p_VJP::Union{Function, DIVJP},
        max_value, min_value, name
) where {T}
    Law{T}(
        GenInputsAndApply(inputs, f),
        GenInputsAndApply(inputs, f_VJP_input),
        GenInputsAndApply(inputs, f_VJP_θ),
        init_cache, callback_freq,
        GenInputsAndApply(inputs, p_VJP),
        CustomVJP(),
        max_value::Real, min_value::Real, name::Symbol
    )
end
function Law{T}( # 4. With precomputed VJP and no on-the-fly VJP computation
        inputs::Union{<: NamedTuple, <: Tuple},
        f::Function, f_VJP_input::Nothing, f_VJP_θ::Nothing,
        init_cache, callback_freq::CB_FREQ_TYPE,
        p_VJP::Union{Function, DIVJP},
        max_value::Real, min_value::Real, name::Symbol
) where {T}
    Law{T}(
        GenInputsAndApply(inputs, f),
        GenInputsAndApply(inputs, emptyVJPWithInputs),
        GenInputsAndApply(inputs, emptyVJPWithInputs),
        init_cache, callback_freq,
        GenInputsAndApply(inputs, p_VJP),
        CustomVJP(),
        max_value, min_value, name
    )
end

### B. Declaration of the law with an affect that directly retrieves the inputs ###
function Law{T}( # 1. No custom VJP
        ::Nothing,
        f::Function, f_VJP_input::Nothing, f_VJP_θ::Nothing,
        init_cache, callback_freq::CB_FREQ_TYPE,
        p_VJP::Nothing,
        max_value::Real, min_value::Real, name::Symbol
) where {T}
    Law{T}(
        f, emptyVJP, emptyVJP,
        init_cache, callback_freq,
        emptyPrepVJP,
        DIVJP(),
        max_value, min_value, name
    )
end
function Law{T}( # 2. With VJP computed on-the-fly
        ::Nothing,
        f::Function, f_VJP_input::Function, f_VJP_θ::Function,
        init_cache, callback_freq::CB_FREQ_TYPE,
        p_VJP::Nothing,
        max_value::Real, min_value::Real, name::Symbol
) where {T}
    Law{T}(
        f, f_VJP_input, f_VJP_θ,
        init_cache, callback_freq,
        emptyPrepVJP,
        CustomVJP(),
        max_value, min_value, name
    )
end
function Law{T}( # 3. With precomputed VJP and possibly some VJP computations on-the-fly
        ::Nothing,
        f::Function, f_VJP_input::Function, f_VJP_θ::Function,
        init_cache, callback_freq::CB_FREQ_TYPE,
        p_VJP::Function,
        max_value::Real, min_value::Real, name::Symbol
) where {T}
    Law{T}(
        f, f_VJP_input, f_VJP_θ,
        init_cache, callback_freq,
        p_VJP,
        CustomVJP(),
        max_value, min_value, name
    )
end
function Law{T}( # 4. With precomputed VJP and no on-the-fly VJP computation
        ::Nothing,
        f::Function, f_VJP_input::Nothing, f_VJP_θ::Nothing,
        init_cache, callback_freq::CB_FREQ_TYPE,
        p_VJP::Function,
        max_value::Real, min_value::Real, name::Symbol
) where {T}
    Law{T}(
        f, emptyVJP, emptyVJP,
        init_cache, callback_freq,
        p_VJP,
        CustomVJP(),
        max_value, min_value, name
    )
end

### C. Wrapper that should be used by the user ###
function Law{T}(; # This method relies on the dispatch above
        inputs = nothing,
        f!, f_VJP_input! = nothing, f_VJP_θ! = nothing,
        init_cache, callback_freq = nothing,
        p_VJP! = nothing,
        max_value = NaN, min_value = NaN, name = :unknown
) where {T}
    Law{T}(
        inputs,
        f!, f_VJP_input!, f_VJP_θ!,
        init_cache, callback_freq,
        p_VJP!,
        max_value, min_value, name
    )
end

function apply_law!(law::Law, cache, simulation, glacier_idx, t, θ)
    law.f(cache, simulation, glacier_idx, t, θ)
end
function law_VJP_input(law::Law, cache, simulation, glacier_idx, t, θ)
    law.f_VJP_input(cache, simulation, glacier_idx, t, θ)
end
function law_VJP_θ(law::Law, cache, simulation, glacier_idx, t, θ)
    law.f_VJP_θ(cache, simulation, glacier_idx, t, θ)
end
function precompute_law_VJP(
        law::Law, cache, vjpsPrepLaw::Nothing, simulation, glacier_idx, t, θ)
    law.p_VJP(cache, vjpsPrepLaw, simulation, glacier_idx, t, θ)
end # Dispatch when p_VJP! is something else than DIVJP()
function precompute_law_VJP(
        law::Law{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any,
            <:GenInputsAndApply{<:Any, DIVJP}, <:Any},
        cache,
        vjpsPrepLaw::AbstractPrepVJP,
        simulation,
        glacier_idx,
        t,
        θ)
    # Dispatch when p_VJP!=DIVJP()
    backend = simulation.parameters.UDE.grad.VJP_method.regressorADBackend
    inputs = generate_inputs(law.p_VJP.inputs, simulation, glacier_idx, t)
    ∂θ = ∂law∂θ!(backend, vjpsPrepLaw, inputs, θ)
    cache.vjp_θ .= ∂θ
    ∂inp = ∂law∂inp!(backend, vjpsPrepLaw, inputs, θ)
    cache.vjp_inp .= ∂inp
end

function init_cache(law::Law, simulation, glacier_idx, θ)
    law.init_cache(simulation, glacier_idx, θ)
end
cache_type(law::Law{CACHE_TYPE}) where {CACHE_TYPE} = CACHE_TYPE

is_callback_law(::Law{<:Any, <:Any, <:Any, <:Any, <:Any, Nothing, <:Any, <:Any}) = false
is_callback_law(::Law{<:Any, <:Any, <:Any, <:Any, <:Any, <:Real, <:Any, <:Any}) = true
function is_precomputable_law_VJP(law::Law{CACHE_TYPE, <:Any, <:Any, <:Any, <:Any, <:Any,
        <:GenInputsAndApply, <:Any}) where {CACHE_TYPE}
    !isa(law.p_VJP.f, typeof(emptyPrepVJPWithInputs))
end
function is_precomputable_law_VJP(law::Law)
    !(isa(law.p_VJP, typeof(emptyPrepVJP)) || isnothing(law.p_VJP))
end

function callback_freq(::Law{<:Any, <:Any, <:Any, <:Any, <:Any, Nothing, <:Any})
    throw("This law does not have callback")
end
function callback_freq(law::Law{<:Any, <:Any, <:Any, <:Any, <:Any, <:Real, <:Any})
    law.callback_freq
end

# Define whether inputs are provided or not through GenInputsAndApply
inputs(law::Law{CACHE_TYPE, <: GenInputsAndApply}) where {CACHE_TYPE} = law.f.inputs
inputs(law::Law) = throw("Inputs are not defined.")
inputs_defined(law::Law{CACHE_TYPE, <: GenInputsAndApply}) where {CACHE_TYPE} = true
inputs_defined(law::Law) = false
apply_law_in_model(law::Law) = !is_callback_law(law)

"""
    ConstantLaw{CACHE_TYPE}(init_cache)

Creates a constant law of type `ConstantLaw{CACHE_TYPE}` that holds a fixed value for the entire simulation.

This is useful to inject glacier-specific or global constants into the simulation without modifying them over time.
The update function is a no-op, and only the `init_cache` function matters.

# Arguments

  - `init_cache::Function`: A function `init_cache(simulation, glacier_idx, θ)::T` that provides the constant value.

# Type Parameters

  - `CACHE_TYPE`: The type of the cache. Must be specified manually and should match the return type of `init_cache`.

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

    function ConstantLaw{CACHE_TYPE}(init_cache) where {CACHE_TYPE}
        return new{CACHE_TYPE, typeof(init_cache)}(init_cache)
    end
end

apply_law!(law::ConstantLaw, cache, simulation, glacier_idx, t, θ) = nothing
law_VJP_input(law::ConstantLaw, cache, simulation, glacier_idx, t, θ) = nothing
law_VJP_θ(law::ConstantLaw, cache, simulation, glacier_idx, t, θ) = nothing
function precompute_law_VJP(law::ConstantLaw, cache, vjpsPrepLaw::AbstractPrepVJP,
        simulation, glacier_idx, t, θ)
    nothing
end
function init_cache(law::ConstantLaw, simulation, glacier_idx, θ)
    law.init_cache(simulation, glacier_idx, θ)
end
cache_type(law::ConstantLaw{CACHE_TYPE}) where {CACHE_TYPE} = CACHE_TYPE
is_callback_law(::ConstantLaw) = false
is_precomputable_law_VJP(law::ConstantLaw) = false
callback_freq(::ConstantLaw) = throw("ConstantLaw doesn't have callback")
inputs(law::ConstantLaw) = (;)
inputs_defined(law::ConstantLaw) = true
apply_law_in_model(law::ConstantLaw) = false

"""
    NullLaw <: AbstractLaw

This struct represents a law that is not used in the iceflow model.
"""
struct NullLaw <: AbstractLaw end

apply_law!(law::NullLaw, cache, simulation, glacier_idx, t, θ) = nothing
law_VJP_input(law::NullLaw, cache, simulation, glacier_idx, t, θ) = nothing
law_VJP_θ(law::NullLaw, cache, simulation, glacier_idx, t, θ) = nothing
function precompute_law_VJP(
        law::NullLaw, cache, vjpsPrepLaw::AbstractPrepVJP, simulation, glacier_idx, t, θ)
    nothing
end
function init_cache(law::NullLaw, simulation, glacier_idx, θ)
    ScalarCacheNoVJP(zeros(Sleipnir.Float))
end
cache_type(law::NullLaw) = ScalarCacheNoVJP
is_callback_law(::NullLaw) = false
is_precomputable_law_VJP(law::NullLaw) = false
callback_freq(::NullLaw) = throw("NullLaw doesn't have callback")
inputs(law::NullLaw) = (;)
inputs_defined(law::NullLaw) = true
apply_law_in_model(law::NullLaw) = false

# Display setup
function repr_cache_type(cType)
    # If cType <: Cache, retrieve the type of cType.value, otherwise use the simple cache type
    if cType <: Cache && hasfield(cType, :value)
        fieldtype(cType, :value)
    else
        cType
    end
end
function repr_eval_law(law::AbstractLaw)
    freq_repr = if apply_law_in_model(law)
        "⟳ "
    elseif callback_freq(law)>0
        "↧$(round(callback_freq(law); digits=3)) "
    else
        "↧@start "
    end
    cType = cache_type(law)
    vjp_repr = if cType <: Cache && hasfield(cType, :vjp_inp) && hasfield(cType, :vjp_θ)
        empty_VJP_input = isa(law.f_VJP_input, typeof(emptyVJP)) ||
                          isa(law.f_VJP_input.f, typeof(emptyVJPWithInputs))
        empty_VJP_θ = isa(law.f_VJP_θ, typeof(emptyVJP)) ||
                      isa(law.f_VJP_θ.f, typeof(emptyVJPWithInputs))
        empty_p_VJP = isa(law.p_VJP, typeof(emptyPrepVJP)) ||
                      isa(law.p_VJP.f, typeof(emptyPrepVJPWithInputs))
        precompDI_repr = (!isa(law.p_VJP, typeof(emptyPrepVJP))) &&
                         isa(law.p_VJP.f, DIVJP) ? " (DI)" : ""
        precomp_repr = is_precomputable_law_VJP(law) ? " ✅ precomputed$(precompDI_repr)" :
                       " ❌ precomputed"
        if empty_VJP_input && empty_VJP_θ && empty_p_VJP
            "auto VJP $precomp_repr"
        else
            "custom VJP $precomp_repr"
        end
    else
        ""
    end
    "($freq_repr $vjp_repr)"
end
function Base.show(io::IO, law::Law{CACHE_TYPE, <: GenInputsAndApply}) where {CACHE_TYPE}
    println(io,
        "$(keys(inputs(law))) -> $(repr_cache_type(cache_type(law)))   $(repr_eval_law(law))")
end
function Base.show(io::IO, law::Law)
    println(io,
        "(simulation, t) -> $(repr_cache_type(cache_type(law)))   $(repr_eval_law(law))")
end
function Base.show(io::IO, law::ConstantLaw)
    println(io, "ConstantLaw -> $(repr_cache_type(cache_type(law)))")
end
Base.show(io::IO, law::NullLaw) = println(io, "NullLaw")
Base.show(io::IO, type::MIME"text/plain", law::AbstractLaw) = Base.show(io, law)
