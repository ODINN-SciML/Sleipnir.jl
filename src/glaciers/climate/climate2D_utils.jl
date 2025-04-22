
###############################################
############  FUNCTIONS   #####################
###############################################

export initialize_glacier_climate!, downscale_2D_climate!, downscale_2D_climate,
    get_cumulative_climate!, get_cumulative_climate, apply_t_cumul_grad!,
    apply_t_grad!, trim_period, partial_year, get_longterm_temps

"""
    initialize_glacier_climate!(glacier::AbstractGlacier, params::Parameters)

Initialize the climate data for a given glacier.

# Arguments
- `glacier::AbstractGlacier`: The glacier object to initialize the climate data for.
- `params::Parameters`: The parameters containing simulation settings and paths.

# Description
This function initializes the climate data for a glacier by:
1. Creating a dummy period based on the simulation time span and step.
2. Loading the raw climate data from a NetCDF file.
3. Calculating the cumulative climate data for the dummy period.
4. Downscaling the cumulative climate data to a 2D grid.
5. Retrieving long-term temperature data for the glacier.
6. Storing the climate data in the glacier object, including raw climate data, cumulative climate data, downscaled 2D climate data, long-term temperatures, average temperatures, and average gradients.
"""
function initialize_glacier_climate!(glacier::AbstractGlacier, params::Parameters)
    dummy_period = partial_year(Day, params.simulation.tspan[1]):Day(1):partial_year(Day, params.simulation.tspan[1] + params.simulation.step)
    raw_climate = RasterStack(joinpath(prepro_dir, params.simulation.rgi_paths[glacier.rgi_id], "raw_climate_$(params.simulation.tspan).nc"))
    if Sleipnir.doublePrec
        raw_climate = convertRasterStackToFloat64(raw_climate)
    end
    climate_step = get_cumulative_climate(raw_climate[At(dummy_period)])
    climate_2D_step = downscale_2D_climate(climate_step, glacier)
    longterm_temps = get_longterm_temps(glacier.rgi_id, params, raw_climate)
    glacier.climate = Climate2D(raw_climate = raw_climate,
                            climate_raw_step = raw_climate[At(dummy_period)],
                            #climate_cum_step = raw_climate.sel(time=dummy_period).sum(),
                            climate_step = climate_step,
                            climate_2D_step = climate_2D_step,
                            longterm_temps = longterm_temps,
                            avg_temps = mean(raw_climate[At(dummy_period)].temp),
                            avg_gradients = mean(raw_climate[At(dummy_period)].gradient))

end

"""
    generate_raw_climate_files(rgi_id::String, simparams::SimulationParameters)

Generate raw climate files for a given RGI (Randolph Glacier Inventory) ID and simulation parameters.

# Arguments
- `rgi_id::String`: The RGI ID for which to generate raw climate files.
- `simparams::SimulationParameters`: The simulation parameters containing the time span and RGI paths.

# Description
This function generates raw climate files for a specified RGI ID if they do not already exist. It retrieves raw climate data, ensures the desired period is covered, crops the data to the desired time period, and saves the raw climate data to disk.

# Details
1. Constructs the path to the RGI directory using the provided `rgi_id` and `simparams`.
2. Checks if the raw climate file for the specified time span already exists.
3. If the file does not exist:
    - Retrieves the raw climate data.
    - Ensures the desired period is covered by the climate data.
    - Crops the climate data to the desired time period.
    - Saves the cropped climate data to disk.
    - Triggers garbage collection to free up memory.
"""
function generate_raw_climate_files(rgi_id::String, simparams::SimulationParameters)
    rgi_path = "" # Initialize RGI path to be accessible outside the try block
    try
        rgi_path = joinpath(prepro_dir, simparams.rgi_paths[rgi_id])
    catch
        @error "RGI path not found for: $rgi_id"
    end

    if !isfile(joinpath(rgi_path, "raw_climate_$(simparams.tspan).nc"))
        println("Getting raw climate data for: ", rgi_id)
        # Get raw climate data for gdir
        tspan_date = partial_year(Day, simparams.tspan[1]):Day(1):partial_year(Day, simparams.tspan[2])
        climate =  get_raw_climate_data(rgi_path)
        # Make sure the desired period is covered by the climate data
        period = trim_period(tspan_date, climate)
        if any((dims(climate, Ti)[begin] <= period[begin]) & any(dims(climate, Ti)[end] >= period[end]))
            climate = climate[At(period)] # Crop desired time period
        else
            @warn "No overlapping period available between climate tspan!"
        end
        # Save raw gdir climate on disk
        write(joinpath(rgi_path, "raw_climate_$(simparams.tspan).nc"), climate)
        GC.gc()
    end
end


"""
    get_cumulative_climate!(climate, period; gradient_bounds=[-0.009, -0.003], default_grad=-0.0065)

Calculate and update the cumulative climate data for a given period.

# Keyword arguments
- `climate::Climate`: The climate object containing raw climate data.
- `period::Period`: The period for which to calculate the cumulative climate data.
- `gradient_bounds::Vector{Float64}`: Optional. The bounds within which to clamp the gradient values. Default is `[-0.009, -0.003]`.
- `default_grad::Float64`: Optional. The default gradient value to use. Default is `-0.0065`.

# Updates
- `climate.climate_raw_step`: The raw climate data for the given period.
- `climate.avg_temps`: The average temperature for the given period.
- `climate.avg_gradients`: The average gradient for the given period.
- `climate.climate_step["prcp"]`: The cumulative precipitation for the given period.
- `climate.climate_step["temp"]`: The cumulative temperature for the given period.
- `climate.climate_step["gradient"]`: The cumulative gradient for the given period.
- `climate.climate_step["avg_temp"]`: The average temperature for the given period.
- `climate.climate_step["avg_gradient"]`: The average gradient for the given period.
- `climate.climate_step["ref_hgt"]`: The reference height from the metadata of the raw climate data.
"""
function get_cumulative_climate!(climate, period, gradient_bounds=[-0.009, -0.003], default_grad=-0.0065)
    climate.climate_raw_step = climate.raw_climate[At(period)]
    climate.avg_temps = mean(climate.climate_raw_step.temp)

    climate.avg_gradients = mean(climate.climate_raw_step.gradient)
    climate.climate_raw_step.temp.data .= max.(climate.climate_raw_step.temp, 0.0) # get PDDs
    climate.climate_raw_step.gradient.data .= clamp.(climate.climate_raw_step.gradient.data, gradient_bounds[1], gradient_bounds[2]) # Clip gradients within plausible values
    climate.climate_step["prcp"] = sum(climate.climate_raw_step.prcp)
    climate.climate_step["temp"] = sum(climate.climate_raw_step.temp)
    climate.climate_step["gradient"] = sum(climate.climate_raw_step.gradient)
    climate.climate_step["avg_temp"] = climate.avg_temps
    climate.climate_step["avg_gradient"] = climate.avg_gradients
    climate.climate_step["ref_hgt"] = metadata(climate.climate_raw_step)["ref_hgt"]
end

"""
    get_cumulative_climate(climate; gradient_bounds=[-0.009, -0.003], default_grad=-0.0065)

Calculate cumulative climate statistics from the given climate data.

# Keyword arguments
- `climate::Climate`: A climate object containing temperature, precipitation, and gradient data.
- `gradient_bounds::Vector{Float64}`: A two-element vector specifying the lower and upper bounds for the gradient values. Defaults to `[-0.009, -0.003]`.
- `default_grad::Float64`: The default gradient value to use. Defaults to `-0.0065`.

# Returns
- `climate_sum::Dict{String, Any}`: A dictionary containing the following keys:
  - `"temp"`: The sum of positive degree days (PDDs) from the temperature data.
  - `"prcp"`: The sum of precipitation data.
  - `"gradient"`: The sum of gradient data, clipped within the specified bounds.
  - `"avg_temp"`: The average temperature.
  - `"avg_gradient"`: The average gradient.
  - `"ref_hgt"`: The reference height from the climate metadata.

# Notes
- The temperature data is modified to only include positive degree-day values (PDDs).
- The gradient data is clipped within the specified bounds to ensure plausible values.
"""
function get_cumulative_climate(climate, gradient_bounds=[-0.009, -0.003], default_grad=-0.0065)
    avg_temp = mean(climate.temp)
    avg_gradient = mean(climate.gradient)
    copy_climate = deepcopy(climate)
    copy_climate.temp.data .= max.(copy_climate.temp.data, 0.0) # get PDDs
    copy_climate.gradient.data .= clamp.(copy_climate.gradient.data, gradient_bounds[1], gradient_bounds[2]) # Clip gradients within plausible values
    climate_sum = Dict("temp" => sum(copy_climate.temp),
                       "prcp" => sum(climate.prcp),
                       "gradient" => sum(copy_climate.gradient),
                       "avg_temp" => avg_temp,
                       "avg_gradient" => avg_gradient,
                       "ref_hgt" => metadata(climate)["ref_hgt"])
    return climate_sum
end


"""
    get_raw_climate_data(rgi_path::String) -> RasterStack

Load raw climate data from a specified path.

# Arguments
- `rgi_path::String`: The file path to the directory containing the climate data file.

# Returns
- `RasterStack`: A `RasterStack` object containing the climate data from the specified file.
"""
function get_raw_climate_data(rgi_path::String)
    climate = RasterStack(joinpath(rgi_path, "climate_historical_daily_W5E5.nc"))
    if Sleipnir.doublePrec
        climate = convertRasterStackToFloat64(climate)
    end
    return climate
end

# TODO: make snow/rain thresholds customizable
"""
    apply_t_cumul_grad!(climate_2D_step::Climate2Dstep, S::Matrix{F}) where {F <: AbstractFloat}

Apply temperature and precipitation gradients based on the positive degree day (PDD) and on the elevation matrix `S` to the climate data in `climate_2D_step`.

# Arguments
- `climate_2D_step::Climate2Dstep`: The climate data structure containing temperature, PDD, gradients, and reference height.
- `S::Matrix{F}`: A matrix of elevations.

# Description
This function updates the temperature and PDD fields in `climate_2D_step` by applying the respective gradients based on the difference between the elevation matrix `S` and the reference height. Negative PDD values are cropped to zero. Additionally, the function adjusts the rain and snow fractions based on the updated temperature values.
"""
function apply_t_cumul_grad!(climate_2D_step::Climate2Dstep, S::Matrix{F}) where {F <: AbstractFloat}
    # We apply the gradients to the temperature
    climate_2D_step.temp .= climate_2D_step.temp .+ climate_2D_step.avg_gradient .* (S .- climate_2D_step.ref_hgt)
    climate_2D_step.PDD .= climate_2D_step.PDD .+ climate_2D_step.gradient .* (S .- climate_2D_step.ref_hgt)
    climate_2D_step.PDD .= ifelse.(climate_2D_step.PDD .< 0.0, 0.0, climate_2D_step.PDD) # Crop negative PDD values

    # We adjust the rain/snow fractions with the updated temperature
    climate_2D_step.snow .= ifelse.(climate_2D_step.temp .> 0.0, 0.0, climate_2D_step.snow)
    climate_2D_step.rain .= ifelse.(climate_2D_step.temp .< 0.0, 0.0, climate_2D_step.rain)
end

"""
    apply_t_grad!(climate::RasterStack, dem::Raster)

Apply temperature gradients to the climate data based on the digital elevation model (DEM).

# Arguments
- `climate::RasterStack`: A `RasterStack` object containing climate data, including temperature and gradient information.
- `dem::Raster`: A `Raster` object representing the digital elevation model (DEM) data.

# Description
This function adjusts the temperature data in the `climate` object by applying the temperature gradients. The adjustment is based on the difference between the mean elevation from the DEM data and a reference height specified in the metadata of the `climate` object.
"""
function apply_t_grad!(climate::RasterStack, dem::Raster)
    # We apply the gradients to the temperature
    climate.temp.data .= climate.temp.data .+ climate.gradient.data .* (mean(dem.data[:]) .- metadata(climate)["ref_hgt"])
end

"""
    downscale_2D_climate!(glacier::Glacier2D)

Update the 2D climate structure for a given glacier by downscaling climate data.

# Arguments
- `glacier::Glacier2D`: The glacier object containing the climate data to be downscaled.

# Description
This function updates the 2D climate structure of the given glacier by:
1. Updating the temperature, PDD (Positive Degree Days), snow, and rain fields in the 2D climate step with the corresponding values from the climate step.
2. Updating the gradients and average gradients in the 2D climate step.
3. Applying temperature gradients and computing the snow/rain fraction for the selected period by reprojecting the current `S` with the `RasterStack` structure.

# Notes
- The function modifies the `glacier` object in place.
"""
function downscale_2D_climate!(glacier::Glacier2D)
    # Update 2D climate structure
    climate = glacier.climate
    climate.climate_2D_step.temp .= climate.climate_step["avg_temp"]
    climate.climate_2D_step.PDD .= climate.climate_step["temp"]
    climate.climate_2D_step.snow .= climate.climate_step["prcp"]
    climate.climate_2D_step.rain .= climate.climate_step["prcp"]
    # Update gradients
    climate.climate_2D_step.gradient = climate.climate_step["gradient"]
    climate.climate_2D_step.avg_gradient = climate.climate_step["avg_gradient"]

    # Apply temperature gradients and compute snow/rain fraction for the selected period
    apply_t_cumul_grad!(climate.climate_2D_step, reshape(glacier.S, size(glacier.S))) # Reproject current S with the RasterStack structure
end

"""
    downscale_2D_climate(climate_step::Dict, glacier::Glacier2D) -> Climate2Dstep

Downscales climate data to a 2D grid based on the provided glacier information.

# Arguments
- `climate_step::Dict`: A dictionary containing climate data for a specific time step. Expected keys are:
  - `"avg_temp"`: Average temperature.
  - `"temp"`: Temperature.
  - `"prcp"`: Precipitation.
  - `"gradient"`: Temperature gradient.
  - `"avg_gradient"`: Average temperature gradient.
  - `"ref_hgt"`: Reference height.
- `glacier::Glacier2D`: A `Glacier2D` object containing glacier data. Expected fields are:
  - `S`: Surface elevation data.
  - `Coords`: A dictionary with keys `"lon"` and `"lat"` for longitude and latitude coordinates.

# Returns
- `Climate2Dstep`: A `Climate2Dstep` object containing the downscaled climate data with fields:
  - `temp`: 2D array of temperature.
  - `PDD`: 2D array of positive degree days.
  - `snow`: 2D array of snow precipitation.
  - `rain`: 2D array of rain precipitation.
  - `gradient`: Temperature gradient.
  - `avg_gradient`: Average temperature gradient.
  - `x`: Longitude coordinates.
  - `y`: Latitude coordinates.
  - `ref_hgt`: Reference height.

# Description
This function creates dummy 2D arrays based on the glacier surface elevation data and applies the climate step data to these arrays. It then constructs a `Climate2Dstep` object with the downscaled climate data and applies temperature gradients to compute the snow/rain fraction for the selected period.
"""
function downscale_2D_climate(climate_step::Dict, glacier::Glacier2D)
    # Create dummy 2D arrays to have a base to apply gradients afterwards
    FT = typeof(glacier.S[1])
    dummy_grid = zeros(size(glacier.S))
    temp_2D = climate_step["avg_temp"] .+ dummy_grid
    PDD_2D = climate_step["temp"] .+ dummy_grid
    snow_2D = climate_step["prcp"] .+ dummy_grid
    rain_2D = climate_step["prcp"] .+ dummy_grid
    climate_2D_step = Climate2Dstep(temp=temp_2D,
                       PDD=PDD_2D,
                       snow=snow_2D,
                       rain=rain_2D,
                       gradient=Float64(climate_step["gradient"]),
                       avg_gradient=Float64(climate_step["avg_gradient"]),
                       x=glacier.Coords["lon"],
                       y=glacier.Coords["lat"],
                       ref_hgt=Float64(climate_step["ref_hgt"]))

    # Apply temperature gradients and compute snow/rain fraction for the selected period
    apply_t_cumul_grad!(climate_2D_step, reshape(glacier.S, size(glacier.S))) # Reproject current S with xarray structure

    return climate_2D_step

end

"""
    downscale_2D_climate(glacier::Glacier2D)

Downscales the climate data for a given 2D glacier.

# Arguments
- `glacier::Glacier2D`: The glacier object containing the climate data to be downscaled.

# Returns
- `climate_2D_step`: The downscaled 2D climate data for the glacier.
"""
function downscale_2D_climate(glacier::Glacier2D)
    climate_2D_step = downscale_2D_climate(glacier.climate.climate_step, glacier)
    return climate_2D_step
end

"""
    trim_period(period, climate)

Adjusts the given `period` to fit within the bounds of the `climate` data, ensuring it aligns with hydrological years.

# Arguments
- `period::UnitRange{Date}`: The initial date range to be trimmed.
- `climate::AbstractArray`: The climate data array, which should have a time dimension `Ti`.

# Returns
- `UnitRange{Date}`: The adjusted date range that fits within the climate data's time bounds.

# Details
- If the start of the climate data is later than the start of the period, the period is adjusted to start from October 1st of the year of the climate data's start.
- If the end of the climate data is earlier than the end of the period, the period is adjusted to end on September 30th of the year of the climate data's end.
"""
function trim_period(period, climate)
    head = dims(climate, Ti)[begin]
    if head > period[begin]
        period = Date(year(head), 10, 1):Day(1):period[end] # make it a hydrological year
    end
    tail = dims(climate, Ti)[end]
    if tail > period[end]
        period = period[1]:Day(1):Date(year(tail), 9, 30) # make it a hydrological year
    end

    return period
end

"""
    partial_year(period::Type{<:Period}, float)

Calculate a partial year date based on a floating-point year value.

# Arguments
- `period::Type{<:Period}`: The type of period to use (e.g., `Month`, `Day`).
- `float::Float64`: The floating-point year value.

# Returns
- `Date`: The calculated date corresponding to the partial year.
"""
function partial_year(period::Type{<:Period}, float)
    _year, Δ = divrem(float, 1)
    year_start = Date(convert(Int,_year))
    year = period((year_start + Year(1)) - year_start)
    partial = period(round(Dates.value(year) * Δ))
    year_start + partial
end

"""
    partial_year(float::Float64) -> Float64

Calculate the partial year value based on the given floating-point number.

# Arguments
- `float::Float64`: A floating-point number representing the fraction of the year.

# Returns
- `Float64`: The calculated partial year value.
"""
partial_year(float) = partial_year(Day, float)

"""
    get_longterm_temps(rgi_id::String, params::Parameters)

Calculate the long-term average temperatures for a given glacier.

# Arguments
- `rgi_id::String`: The RGI (Randolph Glacier Inventory) identifier for the glacier.
- `params::Parameters`: A `Parameters` object containing simulation parameters, including paths to necessary data files.

# Returns
- `longterm_temps`: A vector of long-term average temperatures for the glacier.

# Description
This function reads the gridded data and raw climate data for the specified glacier, applies a temperature gradient correction based on the glacier's topography, and then calculates the long-term average temperatures by grouping the temperature data by year.
"""
function get_longterm_temps(rgi_id::String, params::Parameters)
    rgi_path = joinpath(prepro_dir, params.simulation.rgi_paths[rgi_id])
    glacier_gd = RasterStack(joinpath(rgi_path, "gridded_data.nc"))
    climate = RasterStack(joinpath(rgi_path, "raw_climate_$(params.simulation.tspan).nc"))
    if Sleipnir.doublePrec
        glacier_gd = convertRasterStackToFloat64(glacier_gd)
        climate = convertRasterStackToFloat64(climate)
    end
    apply_t_grad!(climate, glacier_gd.topo)
    longterm_temps = mean.(groupby(climate.temp, Ti=>year)).data
    return longterm_temps
end

"""
    get_longterm_temps(rgi_id::String, params::Parameters, climate::RasterStack) -> Array{Float64}

Calculate the long-term average temperatures for a given glacier.

# Arguments
- `rgi_id::String`: The RGI (Randolph Glacier Inventory) identifier for the glacier.
- `params::Parameters`: A struct containing simulation parameters, including paths to RGI data.
- `climate::RasterStack`: A `RasterStack` object containing climate data.

# Returns
- `Array{Float64}`: An array of long-term average temperatures.

# Description
This function retrieves the gridded data for the specified glacier using its RGI identifier. It then applies a temperature gradient to the climate data based on the glacier's topography. Finally, it calculates the long-term average temperatures by grouping the temperature data by year and computing the mean for each group.
"""
function get_longterm_temps(rgi_id::String, params::Parameters, climate::RasterStack)
    glacier_gd = RasterStack(joinpath(prepro_dir, params.simulation.rgi_paths[rgi_id], "gridded_data.nc"))
    apply_t_grad!(climate, glacier_gd.topo)
    longterm_temps = mean.(groupby(climate.temp, Ti=>year)).data
    return longterm_temps
end


