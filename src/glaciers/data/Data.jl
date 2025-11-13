export tdata

"""
    AbstractData

Abstract type that represents data. Used to implement `ThicknessData` and `SurfaceVelocityData`.
"""
abstract type AbstractData end

include("ThicknessData.jl")
include("SurfaceVelocityMapping.jl")
include("SurfaceVelocityData.jl")

"""
    tdata(data::Nothing)
    tdata(data::ThicknessData)
    tdata(data::Nothing, mapping::MeanDateVelocityMapping)
    tdata(data::SurfaceVelocityData, mapping::MeanDateVelocityMapping)

Retrieve the time steps at which data is available for ice thickness and surface velocity data.
If the provided data is `nothing`, returns an empty vector.
"""
tdata(data::Nothing) = Vector{Sleipnir.Float}() # For non existing data in the fields thicknessData and velocityData of the glaciers
tdata(data::ThicknessData) = isnothing(data.t) ? Vector{Sleipnir.Float}() : data.t
tdata(data::Nothing, mapping::MeanDateVelocityMapping) = Vector{Sleipnir.Float}() # For non existing data in the fields velocityData and velocityData of the glaciers
function tdata(data::SurfaceVelocityData, mapping::MeanDateVelocityMapping)
    isnothing(data.date) ? Vector{Sleipnir.Float}() : datetime_to_floatyear.(data.date)
    # if type==:date
    #     isnothing(data.date) ? Vector{Sleipnir.Float}() : datetime_to_floatyear.(data.date)
    # elseif type==:date1
    #     isnothing(data.date1) ? Vector{Sleipnir.Float}() : datetime_to_floatyear.(data.date1)
    # elseif type==:date2
    #     isnothing(data.date2) ? Vector{Sleipnir.Float}() : datetime_to_floatyear.(data.date2)
    # else
    #     throw("Unknown option type=$(type)")
    # end
end
