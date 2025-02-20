export SurfaceVelocityData

mutable struct SurfaceVelocityData{F <: AbstractFloat} <: AbstractData
    x::Vector{F}
    y::Vector{F}
    lat::Vector{F}
    lon::Vector{F}
    vx::Array{F, 3}
    vy::Array{F, 3}
    vabs
    vx_error
    vy_error
    vabs_error
    date
    date1
    date2
    date_error
    # vx_error::Union{Nothing, Vector{F}} = nothing
    # vy_error::Union{Nothing, Vector{F}} = nothing
    # vabs_error::Union{Nothing, Vector{F}} = nothing
    # date::Vector{Date}
    # date1::Vector{Date}
    # date2::Vector{Date}
    # date_error::Vector{Day}
end

"""


Constructor for ice surface velocity data based on Rabatel et. al (2023). 

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
    vabs::Union{Any, Nothing} = nothing,
    vx_error::Union{Any, Nothing} = nothing,
    vy_error::Union{Any, Nothing} = nothing,
    vabs_error::Union{Any, Nothing} = nothing,
    date::Union{Any, Nothing} = nothing,
    date1::Union{Any, Nothing} = nothing,
    date2::Union{Any, Nothing} = nothing,
    date_error::Union{Any, Nothing} = nothing,
    ) where {F <: AbstractFloat}

    ft = typeof(x[begin])

    return SurfaceVelocityData{ft}(x, y, lat, lon, vx, vy, vabs, vx_error, vy_error, vabs_error, date, date1, date2, date_error)
end


include("surfacevelocitydata_utils.jl")