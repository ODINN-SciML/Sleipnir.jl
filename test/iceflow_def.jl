
# Need these types which are defined in Sleipnir in order to instantiate Results
abstract type IceflowModel <: AbstractModel end
abstract type SIAmodel <: IceflowModel end

# Build a minimalistic iceflow model just to be able to instantiate Results
mutable struct SimpleIceflowModel{R <: Real} <: SIAmodel
    S::Union{Matrix{R}, Nothing} # This is the only required field to build Results
end
