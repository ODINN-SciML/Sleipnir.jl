
export Glacier, Climate

abstract type AbstractGlacier end 

include("Climate.jl")

mutable struct Glacier2D{F <: AbstractFloat, I <: Int} <: AbstractGlacier
    rgi_id::Union{String, Nothing}
    gdir::Union{PyObject, Nothing} 
    climate::Union{Climate, Nothing}
    H₀::Union{Matrix{F}, Nothing}
    S::Union{Matrix{F}, Nothing}
    B::Union{Matrix{F}, Nothing}
    V::Union{Matrix{F}, Nothing}
    slope::Union{Matrix{F}, Nothing}
    dist_border::Union{Matrix{F}, Nothing}
    S_coords::Union{PyObject, Nothing}
    Δx::Union{F, Nothing}
    Δy::Union{F, Nothing}
    nx::Union{I, Nothing}
    ny::Union{I, Nothing}
end

"""
function Glacier(;
    rgi_id::Union{String, Nothing} = nothing,
    gdir::Union{PyObject, Nothing} = nothing,
    climate::Union{Climate, Nothing} = nothing,
    H₀::Union{Matrix{F}, Nothing} = nothing,
    S::Union{Matrix{F}, Nothing} = nothing,
    B::Union{Matrix{F}, Nothing} = nothing,
    V::Union{Matrix{F}, Nothing}= nothing,
    slope::Union{Matrix{F}, Nothing} = nothing,
    dist_border::Union{Matrix{F}, Nothing} = nothing,
    S_coords::Union{PyObject, Nothing} = nothing,
    Δx::Union{F, Nothing} = nothing,
    Δy::Union{F, Nothing} = nothing,
    nx::Union{I, Nothing} = nothing,
    ny::Union{I, Nothing} = nothing
    ) where {F <: AbstractFloat, I <: Int} 

Constructor for empty 2D Glacier object.
"""
function Glacier(;
    rgi_id::Union{String, Nothing} = nothing,
    gdir::Union{PyObject, Nothing} = nothing,
    climate::Union{Climate, Nothing} = nothing,
    H₀::Union{Matrix{F}, Nothing} = nothing,
    S::Union{Matrix{F}, Nothing} = nothing,
    B::Union{Matrix{F}, Nothing} = nothing,
    V::Union{Matrix{F}, Nothing}= nothing,
    slope::Union{Matrix{F}, Nothing} = nothing,
    dist_border::Union{Matrix{F}, Nothing} = nothing,
    S_coords::Union{PyObject, Nothing} = nothing,
    Δx::Union{F, Nothing} = nothing,
    Δy::Union{F, Nothing} = nothing,
    nx::Union{I, Nothing} = nothing,
    ny::Union{I, Nothing} = nothing
    ) where {F <: AbstractFloat, I <: Int} 

    # Define default float and integer type for constructor
    ft = Float64
    it = Int64
    return Glacier2D{ft,it}(rgi_id, gdir, climate, H₀, S, B, V, slope, dist_border, S_coords, Δx, Δy, nx, ny)
end

###############################################
################### UTILS #####################
###############################################

Base.:(==)(a::Glacier2D, b::Glacier2D) = a.rgi_id == b.rgi_id && a.gdir == b.gdir && a.climate == b.climate && 
                                      a.H₀ == b.H₀ && a.S == b.S && a.B == b.B && a.V == b.V &&
                                      a.slope == b.slope && a.dist_border == b.dist_border && a.rgi_id == b.rgi_id &&
                                      a.S_coords == b.S_coords && a.Δx == b.Δx && a.Δy == b.Δy && a.Δx == b.Δx && a.nx == b.nx && a.ny == b.ny

include("glacier_utils.jl")
include("climate_utils.jl")