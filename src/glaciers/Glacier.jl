
export Glacier, Climate

include("Climate.jl")

@kwdef mutable struct Glacier{F <: AbstractFloat, I <: Int} 
    rgi_id::String
    gdir::Union{PyObject, Nothing} 
    climate::Union{Climate, Nothing}
    H₀::Matrix{F}
    S::Matrix{F}
    B::Matrix{F}
    V::Matrix{F}
    slope::Matrix{F}
    dist_border::Matrix{F}
    S_coords::Union{PyObject, Nothing}
    Δx::F
    Δy::F
    nx::I
    ny::I
end

###############################################
################### UTILS #####################
###############################################

Base.:(==)(a::Glacier, b::Glacier) = a.rgi_id == b.rgi_id && a.gdir == b.gdir && a.climate == b.climate && 
                                      a.H₀ == b.H₀ && a.S == b.S && a.B == b.B && a.V == b.V &&
                                      a.slope == b.slope && a.dist_border == b.dist_border && a.rgi_id == b.rgi_id &&
                                      a.S_coords == b.S_coords && a.Δx == b.Δx && a.Δy == b.Δy && a.Δx == b.Δx && a.nx == b.nx && a.ny == b.ny

include("glacier_utils.jl")
include("climate_utils.jl")