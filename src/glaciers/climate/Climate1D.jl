
@kwdef mutable struct Climate1Dstep{F <: AbstractFloat} 
    temp::Vector{F}
    PDD::Vector{F}
    snow::Vector{F}
    rain::Vector{F}
    gradient::Ref{F}
    avg_gradient::Ref{F}
    x::Vector{F}
    y::Vector{F}
    ref_hgt::Ref{F}
end

@kwdef mutable struct Climate1D{F <: AbstractFloat} 
    raw_climate::Py # Raw climate dataset for the whole simulation
    # Buffers to avoid memory allocations
    climate_raw_step::Ref{Py} # Raw climate trimmed for the current step
    climate_step::Ref{Py} # Climate data for the current step
    climate_2D_step::Climate2Dstep # 2D climate data for the current step to feed to the MB model
    longterm_temps::Vector{F} # Longterm temperatures for the ice rheology
    avg_temps::Ref{Py} # Intermediate buffer for computing average temperatures
    avg_gradients::Ref{Py} # Intermediate buffer for computing average gradients
end
