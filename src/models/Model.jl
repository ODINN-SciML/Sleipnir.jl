
export Model, AbstractModel

abstract type AbstractModel end

const AbstractEmptyModel = Union{AbstractModel,Nothing}

# Composite type as a representation of ODINN models
struct Model{IFM <: AbstractEmptyModel, MBM <: AbstractEmptyModel, MLM <: AbstractEmptyModel}
    iceflow::IFM
    mass_balance::MBM
    machine_learning::MLM
end


