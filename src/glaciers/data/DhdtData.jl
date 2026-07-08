export DhdtData

"""
Simple snapshot of mean glacier surface elevation change.
This represents the glacier-wide average dh/dt.
The convention is that negative glacier-wide MB is represented as a negative dh/dt (we assume that the ice density is constant).

Note: if you use an external dataset please make sure that you use the correct unit (surface elevation change vs m.w.e.)
"""
mutable struct DhdtData{F <: AbstractFloat} <: AbstractData
    t::Tuple{F, F}
    dhdt::F
    uncertainty::F

    function DhdtData(t, dhdt, uncertainty = NaN)
        return new{Sleipnir.Float}(t, dhdt, Sleipnir.Float(uncertainty))
    end
end

function Base.:(==)(a::DhdtData, b::DhdtData)
    a.t == b.t && a.dhdt == b.dhdt && a.uncertainty == b.uncertainty
end

function Base.:(≈)(a::DhdtData, b::DhdtData)
    a.t == b.t && safe_approx(a.dhdt, b.dhdt) && safe_approx(a.uncertainty, b.uncertainty)
end
