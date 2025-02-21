export SurfaceVelocityData

mutable struct SurfaceVelocityData{F <: AbstractFloat} <: AbstractData
    x::Union{Vector{F}, Nothing}
    y::Union{Vector{F}, Nothing}
    lat::Union{Vector{F}, Nothing}
    lon::Union{Vector{F}, Nothing}
    vx::Union{Array{F, 3}, Nothing}
    vy::Union{Array{F, 3}, Nothing}
    vabs::Union{Array{F, 3}, Nothing}
    vx_error::Union{Array{F, 1}, Nothing}
    vy_error::Union{Array{F, 1}, Nothing}
    vabs_error::Union{Array{F, 1}, Nothing}
    date::Union{Vector{DateTime}, Nothing}
    date1::Union{Vector{DateTime}, Nothing}
    date2::Union{Vector{DateTime}, Nothing}
    date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing}
end

"""
function SurfaceVelocityData(;
    x::Union{Vector{F}, Nothing} = nothing,
    y::Union{Vector{F}, Nothing} = nothing, 
    lat::Union{Vector{F}, Nothing} = nothing,
    lon::Union{Vector{F}, Nothing} = nothing, 
    vx::Union{Array{F, 3}, Nothing} = nothing,
    vy::Union{Array{F, 3}, Nothing} = nothing,
    vabs::Union{Array{F, 3}, Nothing} = nothing,
    vx_error::Union{Array{F, 1}, Nothing} = nothing,
    vy_error::Union{Array{F, 1}, Nothing} = nothing,
    vabs_error::Union{Array{F, 1}, Nothing} = nothing,
    date::Union{Vector{DateTime}, Nothing} = nothing,
    date1::Union{Vector{DateTime}, Nothing} = nothing,
    date2::Union{Vector{DateTime}, Nothing} = nothing,
    date_error::Union{Vector{Day}, Nothing} = nothing,
    ) where {F <: AbstractFloat}

Constructor for ice surface velocity data based on Rabatel et. al (2023). 


Important remarks:
- Projections in longitude and latitude assume we are working in the north hemisphere. 
  If working with south hemisphere glaciers, this needs to be changed.
- The error in velocity is unique per timestamp, rather than being pixel distributed. 
- The error in the absolute velocities `vabs_error` is overestimated.

References:
    - Rabatel, A., Ducasse, E., Millan, R. & Mouginot, J. 
    Satellite-Derived Annual Glacier Surface Flow Velocity Products for the European Alps, 2015–2021. 
    Data 8, 66 (2023).
"""
function SurfaceVelocityData(;
    x::Union{Vector{F}, Nothing} = nothing,
    y::Union{Vector{F}, Nothing} = nothing, 
    lat::Union{Vector{F}, Nothing} = nothing,
    lon::Union{Vector{F}, Nothing} = nothing, 
    vx::Union{Array{F, 3}, Nothing} = nothing,
    vy::Union{Array{F, 3}, Nothing} = nothing,
    vabs::Union{Array{F, 3}, Nothing} = nothing,
    vx_error::Union{Array{F, 1}, Nothing} = nothing,
    vy_error::Union{Array{F, 1}, Nothing} = nothing,
    vabs_error::Union{Array{F, 1}, Nothing} = nothing,
    date::Union{Vector{DateTime}, Nothing} = nothing,
    date1::Union{Vector{DateTime}, Nothing} = nothing,
    date2::Union{Vector{DateTime}, Nothing} = nothing,
    date_error::Union{Vector{Day}, Vector{Millisecond}, Nothing} = nothing,
    ) where {F <: AbstractFloat}

    ft = typeof(x[begin])

    return SurfaceVelocityData{ft}(x, y, lat, lon, vx, vy, vabs, vx_error, vy_error, vabs_error, date, date1, date2, date_error)
end


include("surfacevelocitydata_utils.jl")