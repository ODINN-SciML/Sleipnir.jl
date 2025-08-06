export ThicknessData

"""
Simple time series of ice thickness data to test transient inversion
"""
mutable struct ThicknessData{F <: AbstractFloat} <: AbstractData
    t::Union{Vector{F}, Nothing}
    H::Union{Vector{Matrix{F}}, Nothing}
end

function ThicknessData(;
    t::Union{Vector{F}, Nothing} = nothing,
    H::Union{Vector{Matrix{F}}, Nothing} = nothing
    ) where {F <: AbstractFloat}
    return ThicknessData{Sleipnir.Float}(t, H)
end

Base.:(==)(a::ThicknessData, b::ThicknessData) = a.t == b.t && a.H == b.H

Base.:(≈)(a::ThicknessData, b::ThicknessData) = safe_approx(a.t, b.t) && safe_approx(a.H, b.H)
