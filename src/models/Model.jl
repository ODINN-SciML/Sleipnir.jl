
export AbstractModel, Model

"""
    AbstractModel

An abstract type that serves as a base for all model types in ODINN.
"""
abstract type AbstractModel end

const AbstractEmptyModel = Union{AbstractModel, Nothing}

"""
    Model{IFM <: AbstractEmptyModel, MBM <: Union{<:AbstractEmptyModel, Vector{<:AbstractEmptyModel}}, TC <: AbstractEmptyModel}

A mutable struct that represents a model with three components: iceflow, mass balance, and machine learning.

    Model(
        iceflow::IFM,
        mass_balance::MBM,
        trainable_components::TC,
    ) where {IFM <: AbstractEmptyModel, MBM <: Union{<:AbstractEmptyModel, Vector{<:AbstractEmptyModel}}, TC <: AbstractEmptyModel}

    Model(;iceflow, mass_balance) = Model(iceflow, mass_balance, nothing)

Initialize Model (no machine learning model).

# Keyword arguments

  - `iceflow::IFM`: Represents the iceflow component, which is an instance of `IFM`.
  - `mass_balance::MBM`: Represents the mass balance component. Either a single `MBM` instance or a `Vector{MBM}` of per-glacier models.
  - `trainable_components::TC`: Represents the trainable components, which is an instance of `TC`.

# Type Parameters

  - `IFM`: A subtype of `AbstractEmptyModel` representing the type of the iceflow model.
  - `MBM`: Either a subtype of `AbstractEmptyModel` (single shared model) or a `Vector` of such subtypes (per-glacier models). The field stores the exact runtime type — no `Union` in the field — which is required for AD compatibility.
  - `TC`: A subtype of `AbstractEmptyModel` representing the type of the trainable components.
"""
mutable struct Model{
    IFM <: AbstractEmptyModel,
    MBM <: Union{<:AbstractEmptyModel, Vector{<:AbstractEmptyModel}},
    TC <: AbstractEmptyModel}
    iceflow::IFM
    mass_balance::MBM
    trainable_components::TC

    function Model(
            iceflow::IFM,
            mass_balance::MBM,
            trainable_components::TC
    ) where {
            IFM <: AbstractEmptyModel,
            MBM <: Union{<:AbstractEmptyModel, Vector{<:AbstractEmptyModel}},
            TC <: AbstractEmptyModel}
        # Parameterize on the concrete runtime types: the `where` bounds above are not
        # guaranteed concrete (e.g. a Union bound has isconcretetype == false), whereas
        # typeof(value) always is, which keeps the fields AD-friendly.
        new{typeof(iceflow), typeof(mass_balance), typeof(trainable_components)}(
            iceflow, mass_balance, trainable_components)
    end
end

function _construct_Model(iceflow::Union{<: AbstractModel, Nothing},
        mass_balance::Union{<: AbstractModel, Nothing}, regressors)
    throw("_construct_Model is not implemented. This is likely because you provided regressors outside of ODINN.")
end
# Since Model is a keyword argument function, we need to define it in Sleipnir
# Otherwise we would have to define it multiple times, making the precompilation impossible with ODINN
function Model(;
        iceflow::Union{<: AbstractModel, Nothing} = nothing,
        mass_balance::Union{<: AbstractModel, Nothing} = nothing,
        regressors::Union{NamedTuple, Nothing} = nothing
)
    if isnothing(regressors)
        Model(iceflow, mass_balance, nothing)
    else
        _construct_Model(iceflow, mass_balance, regressors)
    end
end

"""
    ModelCache{IFC, MBC}

Cache struct that holds the internal state or memory buffers for the components of a `Model`.

Typically used to store per-glacier preallocated buffers or intermediate results
that persist across time steps during simulation.

# Fields

  - `iceflow::IFC`: Cache associated with the iceflow model.
  - `mass_balance::MBC`: Cache associated with the mass balance model.

# Type Parameters

  - `IFC`: Cache type for the iceflow model.
  - `MBC`: Cache type for the mass balance model.
"""
struct ModelCache{IFC, MBC}
    iceflow::IFC
    mass_balance::MBC
end

function init_cache(model::Model, simulation, glacier_idx, θ)
    return ModelCache(
        init_cache(model.iceflow, simulation, glacier_idx, θ),
        # Since mass balance models dont use the "Cache" yet we can just put nothing
        nothing
    )
end
function init_cache(model::Model, simulation, glacier_idx)
    init_cache(model, simulation, glacier_idx, nothing)
end

cache_type(model::Model) = ModelCache{cache_type(model.iceflow), Nothing}

# Display setup
function Base.show(io::IO, type::MIME"text/plain", model::Model)
    println(io, "**** Model ****")
    println(io)
    Base.show(io, type, model.iceflow)
    println(io)
    if model.mass_balance isa AbstractVector
        # Per-glacier mass balance models (e.g. calibrated TImodel1): show a
        # summary plus the first model as a representative sample.
        n = length(model.mass_balance)
        println(io, "Per-glacier mass balance models ($n × $(nameof(eltype(model.mass_balance))))")
        if n > 0
            Base.show(io, type, first(model.mass_balance))
            println(io)
        end
    else
        Base.show(io, type, model.mass_balance)
    end
    println(io)
    if isnothing(model.trainable_components)
        println(io, "No learnable components")
    else
        println(io, "Learnable components")
        Base.show(io, type, model.trainable_components)
        println(io)
    end
    print(io, "***************")
end
