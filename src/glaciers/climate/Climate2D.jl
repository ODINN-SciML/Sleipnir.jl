
export Climate2Dstep, Climate2D

"""
    Climate2Dstep{F <: AbstractFloat}

A mutable struct representing a 2D climate time step with various climate-related parameters.

# Keyword arguments
- `temp::Matrix{F}`: Temperature matrix.
- `PDD::Matrix{F}`: Positive Degree Days matrix.
- `snow::Matrix{F}`: Snowfall matrix.
- `rain::Matrix{F}`: Rainfall matrix.
- `gradient::F`: Gradient value.
- `avg_gradient::F`: Average gradient value.
- `x::Vector{F}`: X-coordinates vector.
- `y::Vector{F}`: Y-coordinates vector.
- `ref_hgt::F`: Reference height.
"""
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

"""
A mutable struct representing a 2D climate for a glacier with various buffers and datasets.

    Climate2D{F <: AbstractFloat}

# Keyword arguments
- `raw_climate::RasterStack`: Raw climate dataset for the whole simulation.
- `climate_raw_step::RasterStack`: Raw climate trimmed for the current step to avoid memory allocations.
- `climate_step::Dict`: Climate data for the current step.
- `climate_2D_step::Climate2Dstep`: 2D climate data for the current step to feed to the mass balance (MB) model.
- `longterm_temps::Vector{F}`: Long-term temperatures for the ice rheology.
- `avg_temps::F`: Intermediate buffer for computing average temperatures.
- `avg_gradients::F`: Intermediate buffer for computing average gradients.
"""
@kwdef mutable struct Climate2D{F <: AbstractFloat}
    raw_climate::RasterStack
    climate_raw_step::RasterStack
    climate_step::Dict
    climate_2D_step::Climate2Dstep
    longterm_temps::Vector{F}
    avg_temps::F
    avg_gradients::F
end

Base.:(==)(a::Climate2D, b::Climate2D) = a.raw_climate == b.raw_climate && a.climate_raw_step == b.climate_raw_step &&
                                      a.climate_step == b.climate_step && a.climate_2D_step == b.climate_2D_step &&
                                      a.longterm_temps == b.longterm_temps && a.avg_temps == b.avg_temps &&
                                      a.avg_gradients == b.avg_gradients

"""
    DummyClimate2D(;
        longterm_temps::Vector{F} = []
    ) where {F <: AbstractFloat}

Dummy climate initialization for very specific use cases where we don't have climate
data and we need to build a minimalistic climate with only a few data.
For the moment it supports only the initialization of the long term temperatures.
It returns a minimalistic Climate2D instance.

Arguments:
- `longterm_temps::Vector{F}`: Long term temperatures.
"""

function DummyClimate2D(;
    longterm_temps::Vector{F} = []
) where {F <: AbstractFloat}
    ras = Raster(rand(X(1:0), Y(1:0), Ti(DateTime(2001):Month(1):DateTime(2002))))
    emptyRasterStack = RasterStack(ras)
    dummyMatrix = [0.;;]
    emptyClimate2Dstep = Climate2Dstep(
        temp = dummyMatrix,
        PDD = dummyMatrix,
        snow = dummyMatrix,
        rain = dummyMatrix,
        gradient = 0.,
        avg_gradient = 0.,
        x = [0.],
        y = [0.],
        ref_hgt = 0.
    )
    return Climate2D(
        raw_climate = emptyRasterStack,
        climate_raw_step = emptyRasterStack,
        climate_step = Dict(),
        climate_2D_step = emptyClimate2Dstep,
        longterm_temps = longterm_temps,
        avg_temps = 0.,
        avg_gradients = 0.,
    )
end
