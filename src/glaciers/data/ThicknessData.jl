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

    ft = Sleipnir.Float

    return ThicknessData{ft}(t, H)
end