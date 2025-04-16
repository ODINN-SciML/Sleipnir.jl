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
- `vx::Union{Vector{Matrix{F}}, Nothing}`: x component of surface velocity.
- `vy::Union{Vector{Matrix{F}}, Nothing}`: y component of surface velocity.
- `vabs::Union{Vector{Matrix{F}}, Nothing}`: Absolute ice surface velocity.
- `vx_error::Union{Vector{F}, Nothing}`: Error in `vx`
- `vy_error::Union{Vector{F}, Nothing}`: Error in `vy`
- `vabs_error::Union{Vector{F}, Nothing}`: Error in `vabs`.
- `date::Union{Vector{DateTime}, Nothing}`: Date of observation (mean of `date1` and `date2`)
- `date1::Union{Vector{DateTime}, Nothing}`: First date of adquisition.
- `date2::Union{Vector{DateTime}, Nothing}`: Second date of adquisition.
- `date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing}`: Error in `date`.
- `glacierGridded::Bool`: Whether the data have been gridded to the glacier grid or not.
"""
mutable struct SurfaceVelocityData{F <: AbstractFloat} <: AbstractData
    x::Union{Vector{F}, Nothing}
    y::Union{Vector{F}, Nothing}
    lat::Union{Vector{F}, Nothing}
    lon::Union{Vector{F}, Nothing}
    vx::Union{Vector{Matrix{F}}, Nothing}
    vy::Union{Vector{Matrix{F}}, Nothing}
    vabs::Union{Vector{Matrix{F}}, Nothing}
    vx_error::Union{Vector{F}, Nothing}
    vy_error::Union{Vector{F}, Nothing}
    vabs_error::Union{Vector{F}, Nothing}
    date::Union{Vector{DateTime}, Nothing}
    date1::Union{Vector{DateTime}, Nothing}
    date2::Union{Vector{DateTime}, Nothing}
    date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing}
    glacierGridded::Bool
end

"""
Constructs `SurfaceVelocityData` using data from Rabatel et. al (2023)
with the given parameters, including default ones.

function SurfaceVelocityData(;
    x::Union{Vector{F}, Nothing} = nothing,
    y::Union{Vector{F}, Nothing} = nothing,
    lat::Union{Vector{F}, Nothing} = nothing,
    lon::Union{Vector{F}, Nothing} = nothing,
    vx::Union{Vector{Matrix{F}}, Nothing} = nothing,
    vy::Union{Vector{Matrix{F}}, Nothing} = nothing,
    vabs::Union{Vector{Matrix{F}}, Nothing} = nothing,
    vx_error::Union{Vector{F}, Nothing} = nothing,
    vy_error::Union{Vector{F}, Nothing} = nothing,
    vabs_error::Union{Vector{F}, Nothing} = nothing,
    date::Union{Vector{DateTime}, Nothing} = nothing,
    date1::Union{Vector{DateTime}, Nothing} = nothing,
    date2::Union{Vector{DateTime}, Nothing} = nothing,
    date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing} = nothing,
    glacierGridded::Bool = false,
) where {F <: AbstractFloat}

Constructor for ice surface velocity data based on Rabatel et. al (2023).


Important remarks:
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
    vx::Union{Vector{Matrix{F}}, Nothing} = nothing,
    vy::Union{Vector{Matrix{F}}, Nothing} = nothing,
    vabs::Union{Vector{Matrix{F}}, Nothing} = nothing,
    vx_error::Union{Vector{F}, Nothing} = nothing,
    vy_error::Union{Vector{F}, Nothing} = nothing,
    vabs_error::Union{Vector{F}, Nothing} = nothing,
    date::Union{Vector{DateTime}, Nothing} = nothing,
    date1::Union{Vector{DateTime}, Nothing} = nothing,
    date2::Union{Vector{DateTime}, Nothing} = nothing,
    date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing} = nothing,
    glacierGridded::Bool = false,
) where {F <: AbstractFloat}

    ft = Sleipnir.Float
    return SurfaceVelocityData{ft}(
        x, y, lat, lon,
        vx, vy, vabs, vx_error, vy_error, vabs_error,
        date, date1, date2, date_error, glacierGridded
    )
end
