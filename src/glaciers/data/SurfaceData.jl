export SurfaceData

"""
    SurfaceData(; t, S, σ)

Dated maps of the surface elevation, for example DEMs, used as observations of the surface.

  - `t::Vector{F}`: Date of each map, in decimal years. The constructor also accepts `Date`s and `DateTime`s.
  - `S::Vector{Matrix{F}}`: Surface elevation (m) of each map on the glacier grid. `NaN` on the
    cells that must not be used.
  - `σ::Vector{Matrix{F}}`: Uncertainty (m, one standard deviation) of each cell given by the
    product, `0` when the product has none and `NaN` where `S` is `NaN`. It does not include
    systematic errors, such as the penetration of the radar into snow and firn.
"""
mutable struct SurfaceData{F <: AbstractFloat} <: AbstractData
    t::Vector{F}
    S::Vector{Matrix{F}}
    σ::Vector{Matrix{F}}
end

function SurfaceData(;
        t::Vector{<:Union{Real, Date, DateTime}},
        S::Vector{Matrix{F}},
        σ::Vector{Matrix{F}}
) where {F <: AbstractFloat}
    return SurfaceData{Sleipnir.Float}(to_floatyear.(t), S, σ)
end

# The maps have `NaN` on the cells that are not used
function Base.:(==)(a::SurfaceData, b::SurfaceData)
    a.t == b.t && isequal(a.S, b.S) && isequal(a.σ, b.σ)
end

function Base.:(≈)(a::SurfaceData, b::SurfaceData)
    a.t ≈ b.t && isapprox(a.S, b.S; nans = true) && isapprox(a.σ, b.σ; nans = true)
end
