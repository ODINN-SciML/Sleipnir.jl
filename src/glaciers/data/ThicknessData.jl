export ThicknessData

"""
Simple time series of ice thickness data to test transient inversion

# Fields

  - `t::Union{Vector{F}, Nothing}`: Date of each observation.
  - `H::Union{Vector{Matrix{F}}, Nothing}`: Observed ice thickness, 0 where it is not observed.
  - `σ::Union{Vector{Matrix{F}}, Nothing}`: Uncertainty of each observed cell, in m, 0 where it
    is not observed. `nothing` if the data has none.
"""
mutable struct ThicknessData{F <: AbstractFloat} <: AbstractData
    t::Union{Vector{F}, Nothing}
    H::Union{Vector{Matrix{F}}, Nothing}
    σ::Union{Vector{Matrix{F}}, Nothing}
end

function ThicknessData(;
        t::Union{Vector{F}, Nothing} = nothing,
        H::Union{Vector{Matrix{F}}, Nothing} = nothing,
        σ::Union{Vector{Matrix{F}}, Nothing} = nothing
) where {F <: AbstractFloat}
    return ThicknessData{Sleipnir.Float}(t, H, σ)
end
function ThicknessData(t::Vector{F}, H::Vector{Matrix{F}}) where {F <: AbstractFloat}
    return ThicknessData{Sleipnir.Float}(t, H, nothing)
end

function Base.:(==)(a::ThicknessData, b::ThicknessData)
    a.t == b.t && a.H == b.H && a.σ == b.σ
end

function Base.:(≈)(a::ThicknessData, b::ThicknessData)
    safe_approx(a.t, b.t) && safe_approx(a.H, b.H) && safe_approx(a.σ, b.σ)
end
