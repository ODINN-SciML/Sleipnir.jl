
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
    ClimateStep{F <: AbstractFloat}

Mutable struct that represents a climate step before downscaling.

# Keyword arguments
- `prcp::F`: Cumulative precipitation for the given period.
- `temp::F`: Cumulative temperature at the reference elevation for the given period.
- `gradient::F`: Cumulative temperature gradient for the given period.
- `avg_temp::F`: Average temperature over the time step.
- `avg_gradient::F`: Average temperature gradient over the time step.
- `ref_hgt::F`: Reference elevation of the raw climate data.
"""
@kwdef mutable struct ClimateStep{F <: AbstractFloat}
    prcp::F
    temp::F
    gradient::F
    avg_temp::F
    avg_gradient::F
    ref_hgt::F
end


Base.:(==)(a::ClimateStep, b::ClimateStep) = a.prcp == b.prcp && a.temp == b.temp &&
                                      a.gradient == b.gradient && a.avg_temp == b.avg_temp &&
                                      a.avg_gradient == b.avg_gradient && a.ref_hgt == b.ref_hgt

"""
A mutable struct representing a 2D climate for a glacier with various buffers and datasets.

    Climate2D{CLIMRAW <: RasterStack, CLIMRAWSTEP <: RasterStack, CLIMSTEP <: ClimateStep, CLIM2DSTEP <: Climate2Dstep, F <: AbstractFloat}

# Fields
- `raw_climate::CLIMRAW`: Raw climate dataset for the whole simulation.
- `climate_raw_step::CLIMRAWSTEP`: Raw climate trimmed for the current step to avoid memory allocations.
- `climate_step::ClimateStep`: Climate data for the current step.
- `climate_2D_step::Climate2Dstep`: 2D climate data for the current step to feed to the mass balance (MB) model.
- `longterm_temps::Vector{F}`: Long-term temperatures for the ice rheology.
- `avg_temps::F`: Intermediate buffer for computing average temperatures.
- `avg_gradients::F`: Intermediate buffer for computing average gradients.
- `ref_hgt::F`: Reference elevation of the raw climate data.

    Climate2D(
        rgi_id,
        params::Parameters,
        S::Matrix{<: AbstractFloat},
        Coords::Dict,
    )

Initialize the climate data given a RGI ID, a matrix of surface elevation and glacier coordinates.

# Arguments
- `rgi_id`: The glacier RGI ID.
- `params::Parameters`: The parameters containing simulation settings and paths.
-  `S::Matrix{<: AbstractFloat}`: Matrix of surface elevation used to initialize the downscaled climate data.
- `Coords::Dict`: Coordinates of the glacier.

# Description
This function initializes the climate data for a glacier by:
1. Creating a dummy period based on the simulation time span and step.
2. Loading the raw climate data from a NetCDF file.
3. Calculating the cumulative climate data for the dummy period.
4. Downscaling the cumulative climate data to a 2D grid.
5. Retrieving long-term temperature data for the glacier.
6. Returning the climate data, including raw climate data, cumulative climate data, downscaled 2D climate data, long-term temperatures, average temperatures, and average gradients.

    Climate2D(
        raw_climate::RasterStack,
        climate_raw_step::RasterStack,
        climate_step::ClimateStep,
        climate_2D_step::Climate2Dstep,
        longterm_temps::Vector{<: AbstractFloat},
        avg_temps::AbstractFloat,
        avg_gradients::AbstractFloat,
        ref_hgt::AbstractFloat,
    )

Initialize the climate data with the fields provided as arguments.
Refer to the list of fields for a complete description of the arguments.
"""
mutable struct Climate2D{CLIMRAW <: RasterStack, CLIMRAWSTEP <: RasterStack, CLIMSTEP <: ClimateStep, CLIM2DSTEP <: Climate2Dstep, F <: AbstractFloat}
    raw_climate::CLIMRAW
    climate_raw_step::CLIMRAWSTEP
    climate_step::CLIMSTEP
    climate_2D_step::CLIM2DSTEP
    longterm_temps::Vector{F}
    avg_temps::F
    avg_gradients::F
    ref_hgt::F

    function Climate2D(
        rgi_id,
        params::Parameters,
        S::Matrix{<: AbstractFloat},
        Coords::Dict,
    )
        dummy_period = partial_year(Day, params.simulation.tspan[1]):Day(1):partial_year(Day, params.simulation.tspan[1] + params.simulation.step)
        raw_climate = RasterStack(joinpath(prepro_dir, params.simulation.rgi_paths[rgi_id], "raw_climate_$(params.simulation.tspan).nc"))
        if Sleipnir.doublePrec
            raw_climate = convertRasterStackToFloat64(raw_climate)
        end
        climate_step = get_cumulative_climate(raw_climate[At(dummy_period)])
        climate_2D_step = downscale_2D_climate(climate_step, S, Coords)
        longterm_temps = get_longterm_temps(rgi_id, params, raw_climate)
        climate_raw_step = raw_climate[At(dummy_period)]
        return new{
                typeof(raw_climate),
                typeof(climate_raw_step),
                typeof(climate_step),
                typeof(climate_2D_step),
                Sleipnir.Float,
            }(
            raw_climate,
            climate_raw_step,
            climate_step,
            climate_2D_step,
            longterm_temps,
            mean(climate_raw_step.temp),
            mean(climate_raw_step.gradient),
            metadata(raw_climate)["ref_hgt"],
        )
    end
    function Climate2D(
        raw_climate::RasterStack,
        climate_raw_step::RasterStack,
        climate_step::ClimateStep,
        climate_2D_step::Climate2Dstep,
        longterm_temps::Vector{<: AbstractFloat},
        avg_temps::AbstractFloat,
        avg_gradients::AbstractFloat,
        ref_hgt::AbstractFloat,
    )
        return new{typeof(raw_climate), typeof(climate_raw_step), typeof(climate_step), typeof(climate_2D_step), Sleipnir.Float}(
            raw_climate,
            climate_raw_step,
            climate_step,
            climate_2D_step,
            longterm_temps,
            avg_temps,
            avg_gradients,
            ref_hgt,
        )
    end
end

Base.:(==)(a::Climate2D, b::Climate2D) = a.raw_climate == b.raw_climate && a.climate_raw_step == b.climate_raw_step &&
                                      a.climate_step == b.climate_step && a.climate_2D_step == b.climate_2D_step &&
                                      a.longterm_temps ≈ b.longterm_temps && a.avg_temps ≈ b.avg_temps &&
                                      a.avg_gradients == b.avg_gradients && a.ref_hgt == b.ref_hgt

diffToDict(a::Climate2D, b::Climate2D) = Dict{Symbol, Bool}(
    :raw_climate => a.raw_climate == b.raw_climate,
    :climate_raw_step => a.climate_raw_step == b.climate_raw_step,
    :climate_step => a.climate_step == b.climate_step,
    :climate_2D_step => a.climate_2D_step == b.climate_2D_step,
    :longterm_temps => a.longterm_temps == b.longterm_temps,
    :avg_temps => a.avg_temps == b.avg_temps,
    :avg_gradients => a.avg_gradients == b.avg_gradients,
    :ref_hgt => a.ref_hgt == b.ref_hgt,
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
    climate_step = ClimateStep(
        gradient = 0.0,
        temp = 0.0,
        avg_temp = 0.0,
        avg_gradient = 0.0,
        ref_hgt = 0.0,
        prcp = 0.0,
    )
    return Climate2D(emptyRasterStack, emptyRasterStack, climate_step, emptyClimate2Dstep, longterm_temps, 0.0, 0.0, 0.0)
end

# TODO: update show with ref_hgt
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
    printstyled(Dates.format(dims(climate.climate_raw_step, Ti)[begin], "yyyy-mm-dd");color=:green)
    print(" → ")
    printstyled(Dates.format(dims(climate.climate_raw_step, Ti)[end], "yyyy-mm-dd");color=:green)
    println(" ]")

    println("  climate_step:")
    print("    gradient = ")
    printstyled(round(climate.climate_step.gradient;digits=3);color=:blue)
    println(" °C/m (sum of clipped gradients)")
    print("    temp = ")
    printstyled(round(climate.climate_step.temp;digits=3);color=:blue)
    println(" °C (sum of positive temperatures)")
    print("    avg_temp = ")
    printstyled(round(climate.climate_step.avg_temp;digits=3);color=:blue)
    println(" °C")
    print("    avg_gradient = ")
    printstyled(round(climate.climate_step.avg_gradient;digits=3);color=:blue)
    println(" °C/m")
    print("    ref_hgt = ")
    printstyled(round(climate.climate_step.ref_hgt;digits=1);color=:blue)
    println(" m")
    print("    prcp = ")
    printstyled(round(climate.climate_step.prcp;digits=1);color=:blue)
    println(" kg/m²")

    print("  longterm_temps: ")
    printstyled("$(typeof(climate.longterm_temps))";color=:yellow)
    print(" with ")
    printstyled("$(length(climate.longterm_temps))";color=:red)
    println(" elements")

    print("  ref_hgt = ")
    printstyled(round(climate.ref_hgt;digits=2);color=:blue)
    println(" m")

    print("  raw_climate: ")
    printstyled("RasterStack";color=:yellow)
    print(" with ")
    printstyled("$(length(climate.raw_climate))";color=:red)
    print(" time steps [ ")
    printstyled(Dates.format(dims(climate.raw_climate, Ti)[begin], "yyyy-mm-dd");color=:green)
    print(" → ")
    printstyled(Dates.format(dims(climate.raw_climate, Ti)[end], "yyyy-mm-dd");color=:green)
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
