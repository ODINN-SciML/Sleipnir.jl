
export Model, AbstractModel

abstract type AbstractModel end

const AbstractEmptyModel = Union{AbstractModel,Nothing}

# Composite type as a representation of ODINN models
mutable struct Model{IFM <: AbstractEmptyModel, MBM <: AbstractEmptyModel, MLM <: AbstractEmptyModel}
    iceflow::Union{IFM, Vector{IFM}}
    mass_balance::Union{MBM, Vector{MBM}}
    machine_learning::MLM
end


