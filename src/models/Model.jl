
export Model, AbstractModel

"""
    AbstractModel

An abstract type that serves as a base for all model types in ODINN.
"""
abstract type AbstractModel end

const AbstractEmptyModel = Union{AbstractModel,Nothing}

"""
A mutable struct that represents a model with three components: iceflow, mass balance, and machine learning.

    Model{IFM <: AbstractEmptyModel, MBM <: AbstractEmptyModel, MLM <: AbstractEmptyModel}

# Keyword arguments
- `iceflow::IFM}`: Represents the iceflow component, which is an instance of `IFM`.
- `mass_balance::Union{MBM, Vector{MBM}}`: Represents the mass balance component, which is an instance of `MBM`.
- `machine_learning::MLM`: Represents the machine learning component, which is an instance of `MLM`.

# Type Parameters
- `IFM`: A subtype of `AbstractEmptyModel` representing the type of the iceflow model.
- `MBM`: A subtype of `AbstractEmptyModel` representing the type of the mass balance model.
- `MLM`: A subtype of `AbstractEmptyModel` representing the type of the machine learning model.
"""
mutable struct Model{IFM <: AbstractEmptyModel, MBM <: AbstractEmptyModel, MLM <: AbstractEmptyModel}
    iceflow::IFM
    mass_balance::MBM
    machine_learning::MLM
end

Model(;iceflow, mass_balance) = Model(iceflow, mass_balance, nothing)

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
        nothing,
    )
end
init_cache(model::Model, simulation, glacier_idx) = init_cache(model, simulation, glacier_idx, nothing)

cache_type(model::Model) = ModelCache{cache_type(model.iceflow), Nothing}
