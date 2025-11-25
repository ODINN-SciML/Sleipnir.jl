
@kwdef mutable struct Climate1Dstep{F <: AbstractFloat}
    temp::Vector{F}
    PDD::Vector{F}
    snow::Vector{F}
    rain::Vector{F}
    gradient::F
    avg_gradient::F
    x::Vector{F}
    y::Vector{F}
    ref_hgt::F
end

function Base.:(==)(a::Climate1Dstep, b::Climate1Dstep)
    a.temp == b.temp && a.PDD == b.PDD &&
        a.snow == b.snow && a.rain == b.rain &&
        a.gradient == b.gradient && a.avg_gradient == b.avg_gradient &&
        a.x == b.x && a.y == b.y && a.ref_hgt == b.ref_hgt
end

@kwdef mutable struct Climate1D{F <: AbstractFloat}
    raw_climate::RasterStack # Raw climate dataset for the whole simulation
    # Buffers to avoid memory allocations
    climate_raw_step::RasterStack # Raw climate trimmed for the current step
    climate_step::Dict # Climate data for the current step
    climate_2D_step::Climate2Dstep # 2D climate data for the current step to feed to the MB model
    longterm_temps::Vector{F} # Longterm temperatures for the ice rheology
    avg_temps::F # Intermediate buffer for computing average temperatures
    avg_gradients::F # Intermediate buffer for computing average gradients
end

function Base.:(==)(a::Climate1D, b::Climate1D)
    a.raw_climate == b.raw_climate && a.climate_raw_step == b.climate_raw_step &&
        a.climate_step == b.climate_step && a.climate_2D_step == b.climate_2D_step &&
        a.longterm_temps == b.longterm_temps && a.avg_temps == b.avg_temps &&
        a.avg_gradients == b.avg_gradients
end
