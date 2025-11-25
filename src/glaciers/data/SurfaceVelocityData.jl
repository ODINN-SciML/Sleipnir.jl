export SurfaceVelocityData

"""
A mutable struct representing a surface velocity data. Notice that all fields can be empty
by providing `nothing` as the default value.

`SurfaceVelocityData{F <: AbstractFloat} <: AbstractData`

# Fields

  - `x::Union{Vector{F}, Nothing}`: Easting of observation.
  - `y::Union{Vector{F}, Nothing}`: Northing of observation.
  - `lat::Union{Vector{F}, Nothing}`: Latitude of observation.
  - `lon::Union{Vector{F}, Nothing}`: Longitude of observation.
  - `vx::Union{Vector{Matrix{F}}, Nothing}`: x / longitudinal component of surface velocity. Positive velocities correspond to the direction of increasing index of the glacier x-coordinate.
  - `vy::Union{Vector{Matrix{F}}, Nothing}`: y / latitudinal component of surface velocity. Positive velocities correspond to the direction of increasing index of the glacier y-coordinate.
  - `vabs::Union{Vector{Matrix{F}}, Nothing}`: Absolute ice surface velocity.
  - `vx_error::Union{Vector{F}, Nothing}`: Error in `vx`
  - `vy_error::Union{Vector{F}, Nothing}`: Error in `vy`
  - `vabs_error::Union{Vector{F}, Nothing}`: Error in `vabs`.
  - `date::Union{Vector{DateTime}, Nothing}`: Date of observation (mean of `date1` and `date2`)
  - `date1::Union{Vector{DateTime}, Nothing}`: First date of acquisition.
  - `date2::Union{Vector{DateTime}, Nothing}`: Second date of acquisition.
  - `date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing}`: Error in `date`.
  - `flag::Union{BitMatrix, Nothing}`: Flag indicating whether a pixel is considered as reliable or not.
  - `isGridGlacierAligned::Bool`: Whether the data have been gridded to the glacier grid or not.
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
    flag::Union{BitMatrix, Nothing}
    isGridGlacierAligned::Bool
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
flag::Union{BitMatrix, Nothing} = nothing,
isGridGlacierAligned::Bool = false,
) where {F <: AbstractFloat}

Constructor for ice surface velocity data based on Rabatel et. al (2023).

Important remarks:

  - Velocities values are reported in m/yr. Positive velocities correspond to the direction of increasing index of the glacier.
    When the glacier is oriented in east-west and south-north (see latitude and coordinate ordering), positive velocities of the ice
    surface velocity correspond to positive east-west and south-north velocity component.
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
        flag::Union{BitMatrix, Nothing} = nothing,
        isGridGlacierAligned::Bool = false
) where {F <: AbstractFloat}
    return SurfaceVelocityData{Sleipnir.Float}(
        x, y, lat, lon,
        vx, vy, vabs, vx_error, vy_error, vabs_error,
        date, date1, date2, date_error, flag, isGridGlacierAligned
    )
end

function Base.:(==)(a::SurfaceVelocityData, b::SurfaceVelocityData)
    a.x == b.x && a.y == b.y && a.lat == b.lat && a.lon == b.lon &&
        a.vx == b.vx && a.vy == b.vy && a.vabs == b.vabs &&
        a.vx_error == b.vx_error && a.vy_error == b.vy_error &&
        a.vabs_error == b.vabs_error &&
        a.date == b.date && a.date1 == b.date1 && a.date2 == b.date2 &&
        a.date_error == b.date_error &&
        a.flag == b.flag &&
        a.isGridGlacierAligned == b.isGridGlacierAligned
end

function Base.:(≈)(a::SurfaceVelocityData, b::SurfaceVelocityData)
    safe_approx(a.x, b.x) && safe_approx(a.y, b.y) && safe_approx(a.lat, b.lat) &&
        safe_approx(a.lon, b.lon) &&
        safe_approx(a.vx, b.vx) && safe_approx(a.vy, b.vy) && safe_approx(a.vabs, b.vabs) &&
        safe_approx(a.vx_error, b.vx_error) && safe_approx(a.vy_error, b.vy_error) &&
        safe_approx(a.vabs_error, b.vabs_error) &&
        safe_approx(a.date, b.date) && safe_approx(a.date1, b.date1) &&
        safe_approx(a.date2, b.date2) && safe_approx(a.date_error, b.date_error) &&
        a.flag == b.flag &&
        a.isGridGlacierAligned == b.isGridGlacierAligned
end
