
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
- `iceflow::Union{IFM, Vector{IFM}}`: Represents the iceflow component, which can be a single instance of `IFM` or a vector of `IFM` instances.
- `mass_balance::Union{MBM, Vector{MBM}}`: Represents the mass balance component, which can be a single instance of `MBM` or a vector of `MBM` instances.
- `machine_learning::MLM`: Represents the machine learning component, which is an instance of `MLM`.

# Type Parameters
- `IFM`: A subtype of `AbstractEmptyModel` representing the type of the iceflow model.
- `MBM`: A subtype of `AbstractEmptyModel` representing the type of the mass balance model.
- `MLM`: A subtype of `AbstractEmptyModel` representing the type of the machine learning model.
"""
mutable struct Model{IFM <: AbstractEmptyModel, MBM <: AbstractEmptyModel, MLM <: AbstractEmptyModel}
    iceflow::Union{IFM, Vector{IFM}}
    mass_balance::Union{MBM, Vector{MBM}}
    machine_learning::MLM
end

Model(;iceflow, mass_balance) = Model(iceflow, mass_balance, nothing)

struct ModelCache{IFC, MBC}
    iceflow::IFC
    mass_balance::MBC
end

function init_cache(model::Model, simulation, glacier_idx, θ)
    return ModelCache(
        init_cache(model.iceflow, simulation, glacier_idx, θ),
        nothing, # Since mass balance models dont use the "Cache" yet i just put nothing
    )
end

cache_type(model::Model) = ModelCache{cache_type(model.iceflow), Nothing}
