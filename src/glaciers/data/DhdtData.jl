export DhdtData

"""
Simple snapshot of mean glacier surface elevation change.
Negative MB is represented as a negative dhdt.

Note: if you use an external dataset please make sure that you use the correct unit (surface elevation change vs m.w.e.)
"""
mutable struct DhdtData{F <: AbstractFloat} <: AbstractData
    t::Tuple{F, F}
    dhdt::F

    function DhdtData(t, H)
        return new{Sleipnir.Float}(t, H)
    end
end

Base.:(==)(a::DhdtData, b::DhdtData) = a.t == b.t && a.dhdt == b.dhdt

function Base.:(≈)(a::DhdtData, b::DhdtData)
    safe_approx(a.t, b.t) && safe_approx(a.dhdt, b.dhdt)
end
