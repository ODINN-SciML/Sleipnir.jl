
export Model

abstract type AbstractModel end

# Composite type as a representation of ODINN models
struct Model{M <: AbstractModel}
    iceflow::M
    mass_balance::M
    machine_learning::M
end


