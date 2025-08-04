
export Climate2Dstep, Climate2D

"""
    Climate2Draw{F <: AbstractFloat}

An immutable struct representing raw climate data loaded from OGGM.
This struct represents the data contained in the netCDF file and it ensures type
stability as the concrete type of a `RasterStack` depends on the loaded data.

# Keyword arguments
- `prcp::Vector{F}`: Vector of the total daily precipitation amount.
- `temp::Vector{F}`: Vector of 2m daily temperature at the reference elevation.
- `gradient::Vector{F}`: Vector of daily temperature gradient.
- `ref_hgt::F`: Reference elevation of the raw climate data.
- `Ti::Vector{DateTime}`: Vector of `DateTime` representing the time step of each of the raw climate data.
"""
struct Climate2Draw{F <: AbstractFloat}
    prcp::Vector{F}
    temp::Vector{F}
    gradient::Vector{F}
    ref_hgt::F
    Ti::Vector{DateTime}

    function Climate2Draw(raw_climate_rasterstack)
        return Climate2Draw(
            replace(raw_climate_rasterstack.prcp.data, missing => NaN),
            replace(raw_climate_rasterstack.temp.data, missing => NaN),
            replace(raw_climate_rasterstack.gradient.data, missing => NaN),
            Sleipnir.Float(metadata(raw_climate_rasterstack)["ref_hgt"]),
            dims(raw_climate_rasterstack, Ti).val.data,
        )
    end
    function Climate2Draw(
        prcp::Vector{F} = Vector{Sleipnir.Float}([]),
        temp::Vector{F} = Vector{Sleipnir.Float}([]),
        gradient::Vector{F} = Vector{Sleipnir.Float}([]),
        ref_hgt::F = 0.0,
        Ti::Vector{DateTime} = Vector{DateTime}([]),
    ) where {F <: AbstractFloat}
        N = length(Ti)
        @assert N==length(prcp) "Number of elements in prcp should match number of time steps Ti."
        @assert N==length(temp) "Number of elements in temp should match number of time steps Ti."
        @assert N==length(gradient) "Number of elements in gradient should match number of time steps Ti."
        new{Sleipnir.Float}(
            prcp, temp, gradient, ref_hgt, Ti
        )
    end
end

Base.:(==)(a::Climate2Draw, b::Climate2Draw) = a.prcp == b.prcp && a.temp == b.temp &&
                                      a.gradient == b.gradient && a.ref_hgt == b.ref_hgt &&
                                      a.Ti == b.Ti

Base.length(a::Climate2Draw) = length(a.Ti)

diffToDict(a::Climate2Draw, b::Climate2Draw) = Dict{Symbol, Bool}(
    :prcp => a.prcp == b.prcp,
    :temp => a.temp == b.temp,
    :gradient => a.gradient == b.gradient,
    :ref_hgt => a.ref_hgt == b.ref_hgt,
    :Ti => a.Ti == b.Ti,
)

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
- `raw_climate::Climate2Draw`: Raw climate dataset for the whole simulation.
- `climate_raw_step::Climate2Draw`: Raw climate trimmed for the current step to avoid memory allocations.
- `climate_step::Dict`: Climate data for the current step.
- `climate_2D_step::Climate2Dstep`: 2D climate data for the current step to feed to the mass balance (MB) model.
- `longterm_temps::Vector{F}`: Long-term temperatures for the ice rheology.
- `avg_temps::F`: Intermediate buffer for computing average temperatures.
- `avg_gradients::F`: Intermediate buffer for computing average gradients.
"""
@kwdef mutable struct Climate2D{CLIMRAW <: Climate2Draw, F <: AbstractFloat}
    raw_climate::CLIMRAW
    climate_raw_step::CLIMRAW
    climate_step::Dict
    climate_2D_step::Climate2Dstep
    longterm_temps::Vector{F}
    avg_temps::F
    avg_gradients::F
end

Base.:(==)(a::Climate2D, b::Climate2D) = a.raw_climate == b.raw_climate && a.climate_raw_step == b.climate_raw_step &&
                                      a.climate_step == b.climate_step && a.climate_2D_step == b.climate_2D_step &&
                                      a.longterm_temps ≈ b.longterm_temps && a.avg_temps ≈ b.avg_temps &&
                                      a.avg_gradients == b.avg_gradients

diffToDict(a::Climate2D, b::Climate2D) = Dict{Symbol, Bool}(
    :raw_climate => a.raw_climate == b.raw_climate,
    :climate_raw_step => a.climate_raw_step == b.climate_raw_step,
    :climate_step => a.climate_step == b.climate_step,
    :climate_2D_step => a.climate_2D_step == b.climate_2D_step,
    :longterm_temps => a.longterm_temps == b.longterm_temps,
    :avg_temps => a.avg_temps == b.avg_temps,
    :avg_gradients => a.avg_gradients == b.avg_gradients,
)

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
    longterm_temps::Vector{F} = Vector{Sleipnir.Float}([])
) where {F <: AbstractFloat}
    emptyClimate2Draw = Climate2Draw()
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
    return Climate2D{typeof(emptyClimate2Draw), Sleipnir.Float}(
        raw_climate = emptyClimate2Draw,
        climate_raw_step = emptyClimate2Draw,
        climate_step = Dict(
            "gradient" => 0.0,
            "temp" => 0.0,
            "avg_temp" => 0.0,
            "avg_gradient" => 0.0,
            "ref_hgt" => 0.0,
            "prcp" => 0.0
        ),
        climate_2D_step = emptyClimate2Dstep,
        longterm_temps = longterm_temps,
        avg_temps = 0.,
        avg_gradients = 0.,
    )
end

# Display setup
function Base.show(io::IO, type::MIME"text/plain", climate::Climate2D)
    Base.show(io, climate)
end
function Base.show(io::IO, climate::Climate2D)
    printstyled("Climate2D\n";color=:yellow)

    print("  avg_gradients = ")
    printstyled(round(climate.avg_gradients; digits=4); color=:blue)
    println(" °C/m")

    print("  avg_temps = ")
    printstyled(round(climate.avg_temps; digits=2);color=:blue)
    println(" °C")

    print("  climate_2D_step: ")
    printstyled("Climate2Dstep";color=:yellow)
    print(" with a ")
    printstyled(size(climate.climate_2D_step.temp); color=:red)
    println(" grid")

    print("  climate_raw_step: ")
    printstyled("RasterStack";color=:yellow)
    print(" with ")
    printstyled(length(climate.climate_raw_step);color=:red)
    print(" time steps [ ")
    printstyled(Dates.format(climate.climate_raw_step.Ti[begin], "yyyy-mm-dd");color=:green)
    print(" → ")
    printstyled(Dates.format(climate.climate_raw_step.Ti[end], "yyyy-mm-dd");color=:green)
    println(" ]")

    println("  climate_step:")
    print("    gradient = ")
    printstyled(round(climate.climate_step["gradient"];digits=3);color=:blue)
    println(" °C/m (sum of clipped gradients)")
    print("    temp = ")
    printstyled(round(climate.climate_step["temp"];digits=3);color=:blue)
    println(" °C (sum of positive temperatures)")
    print("    avg_temp = ")
    printstyled(round(climate.climate_step["avg_temp"];digits=3);color=:blue)
    println(" °C")
    print("    avg_gradient = ")
    printstyled(round(climate.climate_step["avg_gradient"];digits=3);color=:blue)
    println(" °C/m")
    print("    ref_hgt = ")
    printstyled(round(climate.climate_step["ref_hgt"];digits=1);color=:blue)
    println(" m")
    print("    prcp = ")
    printstyled(round(climate.climate_step["prcp"];digits=1);color=:blue)
    println(" kg/m²")

    print("  longterm_temps: ")
    printstyled("$(typeof(climate.longterm_temps))";color=:yellow)
    print(" with ")
    printstyled("$(length(climate.longterm_temps))";color=:red)
    println(" elements")

    print("  raw_climate: ")
    printstyled("RasterStack";color=:yellow)
    print(" with ")
    printstyled("$(length(climate.raw_climate))";color=:red)
    print(" time steps [ ")
    printstyled(Dates.format(climate.raw_climate.Ti[begin], "yyyy-mm-dd");color=:green)
    print(" → ")
    printstyled(Dates.format(climate.raw_climate.Ti[end], "yyyy-mm-dd");color=:green)
    print(" ]")
end

function Base.show(io::IO, type::MIME"text/plain", climate_step::Climate2Dstep)
    Base.show(io, climate_step)
end
function Base.show(io::IO, climate_step::Climate2Dstep)
    printstyled("Climate2Dstep";color=:yellow)
    print(" with a ")
    printstyled(size(climate_step.temp); color=:red)
    println(" grid")

    print("  avg_gradient = ")
    printstyled(round(climate_step.avg_gradient; digits=4); color=:blue)
    println(" °C/m")

    print("  gradient = ")
    printstyled(round(climate_step.gradient; digits=4); color=:blue)
    println(" °C/m (sum of clipped gradients)")

    print("  ref_hgt = ")
    printstyled(round(climate_step.ref_hgt; digits=1); color=:blue)
    println(" m")

    print("  x: ")
    printstyled(typeof(climate_step.x);color=:yellow)
    print(" in [ ")
    printstyled("$(round(minimum(climate_step.x);digits=6))° → $(round(maximum(climate_step.x);digits=6))°"; color=:blue)
    println(" ]")

    print("  y: ")
    printstyled(typeof(climate_step.y);color=:yellow)
    print(" in [ ")
    printstyled("$(round(minimum(climate_step.y);digits=6))° → $(round(maximum(climate_step.y);digits=6))°"; color=:blue)
    println(" ]")

    println("  min,mean,max:")

    print("    PDD: ")
    printstyled("$(round(minimum(climate_step.PDD); digits=1)) $(round(mean(climate_step.PDD); digits=1)) $(round(maximum(climate_step.PDD); digits=1))\n"; color=:blue)

    print("    rain: ")
    printstyled("$(round(minimum(climate_step.rain); digits=1)) $(round(mean(climate_step.rain); digits=1)) $(round(maximum(climate_step.rain); digits=1))\n"; color=:blue)

    print("    snow: ")
    printstyled("$(round(minimum(climate_step.snow); digits=1)) $(round(mean(climate_step.snow); digits=1)) $(round(maximum(climate_step.snow); digits=1))\n"; color=:blue)

    print("    temp: ")
    printstyled("$(round(minimum(climate_step.temp); digits=1)) $(round(mean(climate_step.temp); digits=1)) $(round(maximum(climate_step.temp); digits=1))\n"; color=:blue)
end
