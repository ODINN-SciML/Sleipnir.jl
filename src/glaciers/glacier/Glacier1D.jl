
export Glacier1D, Climate1D, AbstractGlacier

abstract type AbstractGlacier end

include("../climate/Climate1D.jl")

mutable struct Glacier1D{F <: AbstractFloat, I <: Integer} <: AbstractGlacier
    rgi_id::String
    climate::Union{Climate1D, Nothing}
    H₀::Vector{F}
    S::Vector{F}
    B::Vector{F}
    V::Vector{F}
    A::Union{F, Nothing}
    C::Union{F, Nothing}
    n::Union{F, Nothing}
    w₀::Union{Vector{F}, Nothing}
    λ::Union{Vector{F}, Nothing}
    slope::Vector{F}
    dist_border::Vector{F}
    Coords::Dict{String, Vector{Float64}}
    Δx::F
    Δy::F
    nx::I
    ny::I
end

"""
function Glacier1D(;
    rgi_id::String = "",
    climate::Union{Climate1D, Nothing} = nothing,
    H₀::Vector{F} = Vector{Sleipnir.Float}([]),
    S::Vector{F} = Vector{Sleipnir.Float}([]),
    B::Vector{F} = Vector{Sleipnir.Float}([]),
    V::Vector{F}= Vector{Sleipnir.Float}([]),
    A::Union{F, Nothing} = nothing,
    C::Union{F, Nothing} = nothing,
    n::Union{F, Nothing} = nothing,
    w₀::Union{Vector{F}, Nothing} = nothing,
    λ::Union{Vector{F}, Nothing} = nothing,
    slope::Vector{F} = Vector{Sleipnir.Float}([]),
    dist_border::Vector{F} = Vector{Sleipnir.Float}([]),
    Coords::Dict{String, Vector{Float64}} = Dict{String, Vector{Float64}}("lon" => [], "lat" => []),
    Δx::F = 0,
    Δy::F = 0,
    nx::I = 0,
    ny::I = 0,
    ) where {F <: AbstractFloat, I <: Integer}

Constructor for empty 2D Glacier object.
"""
function Glacier1D(;
    rgi_id::String = "",
    climate::Union{Climate1D, Nothing} = nothing,
    H₀::Vector{F} = Vector{Sleipnir.Float}([]),
    S::Vector{F} = Vector{Sleipnir.Float}([]),
    B::Vector{F} = Vector{Sleipnir.Float}([]),
    V::Vector{F} = Vector{Sleipnir.Float}([]),
    A::Union{F, Nothing} = nothing,
    C::Union{F, Nothing} = nothing,
    n::Union{F, Nothing} = nothing,
    w₀::Union{Vector{F}, Nothing} = nothing,
    λ::Union{Vector{F}, Nothing} = nothing,
    slope::Vector{F} = Vector{Sleipnir.Float}([]),
    dist_border::Vector{F} = Vector{Sleipnir.Float}([]),
    Coords::Dict{String, Vector{Float64}} = Dict{String, Vector{Float64}}("lon" => [], "lat" => []),
    Δx::F = 0,
    Δy::F = 0,
    nx::I = 0,
    ny::I = 0,
    ) where {F <: AbstractFloat, I <: Integer}

    # Define default float and integer type for constructor
    ft = Sleipnir.Float
    it = Sleipnir.Int

    return Glacier1D{ft,it}(rgi_id, gdir, climate, H₀, S, B, V, A, C, n, w₀, λ, slope, dist_border, Coords, Δx, Δy, nx, ny)
end

###############################################
################### UTILS #####################
###############################################

Base.:(==)(a::Glacier1D, b::Glacier1D) = a.rgi_id == b.rgi_id && a.gdir == b.gdir && a.climate == b.climate &&
                                      a.H₀ == b.H₀ && a.S == b.S && a.B == b.B && a.V == b.V &&
                                      a.A == b.A && a.C == b.C && a.n == b.n && a.w₀ == b.w₀ && a.λ == b.λ &&
                                      a.slope == b.slope && a.dist_border == b.dist_border && a.rgi_id == b.rgi_id &&
                                      a.Coords == b.Coords && a.Δx == b.Δx && a.Δy == b.Δy && a.Δx == b.Δx && a.nx == b.nx && a.ny == b.ny

include("glacier1D_utils.jl")
include("../climate/climate1D_utils.jl")