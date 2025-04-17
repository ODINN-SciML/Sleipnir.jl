export AbstractData

"""
    AbstractData

Abstract type that represents data. Used to implement `ThicknessData` and `SurfaceVelocityData`.
"""
abstract type AbstractData end

include("ThicknessData.jl")
include("SurfaceVelocityMapping.jl")
include("SurfaceVelocityData.jl")
