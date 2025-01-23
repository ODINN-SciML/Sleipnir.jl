
export Climate2Dstep, Climate2D

@kwdef mutable struct Climate2Dstep{F <: AbstractFloat}
    temp::Matrix{F}
    PDD::Matrix{F}
    snow::Matrix{F}
    rain::Matrix{F}
    gradient::F
    avg_gradient::F
    x::Vector{F}
    y::Vector{F}
    ref_hgt::F
end

Base.:(==)(a::Climate2Dstep, b::Climate2Dstep) = a.temp == b.temp && a.PDD == b.PDD &&
                                      a.snow == b.snow && a.rain == b.rain &&
                                      a.gradient == b.gradient && a.avg_gradient == b.avg_gradient &&
                                      a.x == b.x && a.y == b.y && a.ref_hgt == b.ref_hgt

@kwdef mutable struct Climate2D{F <: AbstractFloat}
    raw_climate::RasterStack # Raw climate dataset for the whole simulation
    # Buffers to avoid memory allocations
    climate_raw_step::RasterStack # Raw climate trimmed for the current step
    climate_step::Dict # Climate data for the current step
    climate_2D_step::Climate2Dstep # 2D climate data for the current step to feed to the MB model
    longterm_temps::Vector{F} # Longterm temperatures for the ice rheology
    avg_temps::F # Intermediate buffer for computing average temperatures
    avg_gradients::F # Intermediate buffer for computing average gradients
end

Base.:(==)(a::Climate2D, b::Climate2D) = a.raw_climate == b.raw_climate && a.climate_raw_step == b.climate_raw_step &&
                                      a.climate_step == b.climate_step && a.climate_2D_step == b.climate_2D_step &&
                                      a.longterm_temps == b.longterm_temps && a.avg_temps == b.avg_temps &&
                                      a.avg_gradients == b.avg_gradients
