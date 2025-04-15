export SurfaceVelocityData

include("surfacevelocitydata_utils.jl")

"""
A mutable struct representing a surface velocity data. Notice that all fields can be empty
by providing `nothing` as the default value.

`SurfaceVelocityData{F <: AbstractFloat} <: AbstractData`

# Fields
 - `x::Union{Vector{F}, Nothing}`: Easting of observation.
 - `y::Union{Vector{F}, Nothing}`: Northing of observation.
 - `lat::Union{Vector{F}, Nothing}`: Latitude of observation.
 - `lon::Union{Vector{F}, Nothing}`: Longitude of observation.
 - `vx::Union{Array{Union{Missing, F}, 3}, Nothing}`: x component of surface velocity.
 - `vy::Union{Array{Union{Missing, F}, 3}, Nothing}`: y component of surface velocity.
 - `vabs::Union{Array{Union{Missing, F}, 3}, Nothing}`: Absolute ice surface velocity.
 - `vx_error::Union{Array{F, 1}, Nothing}`: Error in `vx`
 - `vy_error::Union{Array{F, 1}, Nothing}`: Error in `vy`
 - `vabs_error::Union{Array{F, 1}, Nothing}`: Error in `vabs`.
 - `date::Union{Vector{DateTime}, Nothing}`: Date of observation (mean of `date1` and `date2`)
 - `date1::Union{Vector{DateTime}, Nothing}`: First date of adquisition.
 - `date2::Union{Vector{DateTime}, Nothing}`: Second date of adquisition.
 - `date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing}`: Error in `date`.
"""
mutable struct SurfaceVelocityData{F <: AbstractFloat} <: AbstractData
    x::Union{Vector{F}, Nothing}
    y::Union{Vector{F}, Nothing}
    lat::Union{Vector{F}, Nothing}
    lon::Union{Vector{F}, Nothing}
    vx::Union{Array{Union{Missing, F}, 3}, Nothing}
    vy::Union{Array{Union{Missing, F}, 3}, Nothing}
    vabs::Union{Array{Union{Missing, F}, 3}, Nothing}
    vx_error::Union{Array{F, 1}, Nothing}
    vy_error::Union{Array{F, 1}, Nothing}
    vabs_error::Union{Array{F, 1}, Nothing}
    date::Union{Vector{DateTime}, Nothing}
    date1::Union{Vector{DateTime}, Nothing}
    date2::Union{Vector{DateTime}, Nothing}
    date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing}
end

"""
Constructs `SurfaceVelocityData` using data from Rabatel et. al (2023)
with the given parameters, including default ones.

function SurfaceVelocityData(;
    x::Union{Vector{F}, Nothing} = nothing,
    y::Union{Vector{F}, Nothing} = nothing,
    lat::Union{Vector{F}, Nothing} = nothing,
    lon::Union{Vector{F}, Nothing} = nothing,
    vx::Union{Array{Union{Missing, F}, 3}, Nothing} = nothing,
    vy::Union{Array{Union{Missing, F}, 3}, Nothing} = nothing,
    vabs::Union{Array{Union{Missing, F}, 3}, Nothing} = nothing,
    vx_error::Union{Array{F, 1}, Nothing} = nothing,
    vy_error::Union{Array{F, 1}, Nothing} = nothing,
    vabs_error::Union{Array{F, 1}, Nothing} = nothing,
    date::Union{Vector{DateTime}, Nothing} = nothing,
    date1::Union{Vector{DateTime}, Nothing} = nothing,
    date2::Union{Vector{DateTime}, Nothing} = nothing,
    date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing} = nothing,
    ) where {F <: AbstractFloat}

Constructor for ice surface velocity data based on Rabatel et. al (2023).


Important remarks:
- Projections in longitude and latitude assume we are working in the north hemisphere.
  If working with south hemisphere glaciers, this needs to be changed.
- The error in velocity is unique per timestamp, rather than being pixel distributed.
- The error in the absolute velocities `vabs_error` is overestimated.

References:
    - Rabatel, A., Ducasse, E., Millan, R. & Mouginot, J.
    Satellite-Derived Annual Glacier Surface Flow Velocity Products for the European Alps,
    2015–2021.
    Data 8, 66 (2023).
"""
function SurfaceVelocityData(;
    x::Union{Vector{F}, Nothing} = nothing,
    y::Union{Vector{F}, Nothing} = nothing,
    lat::Union{Vector{F}, Nothing} = nothing,
    lon::Union{Vector{F}, Nothing} = nothing,
    vx::Union{Array{Union{Missing, F}, 3}, Nothing} = nothing,
    vy::Union{Array{Union{Missing, F}, 3}, Nothing} = nothing,
    vabs::Union{Array{Union{Missing, F}, 3}, Nothing} = nothing,
    vx_error::Union{Array{F, 1}, Nothing} = nothing,
    vy_error::Union{Array{F, 1}, Nothing} = nothing,
    vabs_error::Union{Array{F, 1}, Nothing} = nothing,
    date::Union{Vector{DateTime}, Nothing} = nothing,
    date1::Union{Vector{DateTime}, Nothing} = nothing,
    date2::Union{Vector{DateTime}, Nothing} = nothing,
    date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing} = nothing,
    ) where {F <: AbstractFloat}

    ft = Sleipnir.Float

    return SurfaceVelocityData{ft}(
        x, y, lat, lon,
        vx, vy, vabs, vx_error, vy_error, vabs_error,
        date, date1, date2, date_error
    )
end
