
###############################################
############  FUNCTIONS   #####################
###############################################

export downscale_2D_climate!, downscale_2D_climate,
       get_cumulative_climate!, get_cumulative_climate,
       apply_t_grad!, trim_period, partial_year, get_longterm_temps,
       get_winter_prcp_factor

function _aggregate_raw_layer(climate_raw_step::RasterStack, layer::Symbol; reducer = sum)
    if hasproperty(climate_raw_step, layer)
        return Sleipnir.Float(round(
            reducer(getproperty(climate_raw_step, layer)); digits = 8))
    end
    return Sleipnir.Float(0.0)
end

# Effective period (decimal years) the raw climate series must cover. When the mass
# balance model is calibrated against geodetic observations, the series must also span
# the Hugonnet calibration window, not just the simulation tspan.
function _climate_period(simparams::SimulationParameters)
    tspan = simparams.tspan
    if simparams.use_MB && simparams.calibrate_MB
        return (min(tspan[1], HUGONNET_CLIMATE_PERIOD[1]),
            max(tspan[2], HUGONNET_CLIMATE_PERIOD[2]))
    end
    return tspan
end

function _slice_climate_between_dates(
        climate::RasterStack, start_date::Date, end_date::Date)
    time_axis = collect(dims(climate, Ti))
    selected = filter(t -> start_date <= Date(t) <= end_date, time_axis)
    isempty(selected) &&
        throw(ArgumentError("No climate timesteps found between $(start_date) and $(end_date)."))
    return climate[At(selected)]
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

      + Retrieves the raw climate data.
      + Ensures the desired period is covered by the climate data.
      + Crops the climate data to the desired time period.
      + Saves the cropped climate data to disk. # Initialize RGI path to be accessible outside the try block
      + Triggers garbage collection to free up memory.
"""
function generate_raw_climate_files(rgi_id::String, simparams::SimulationParameters)
    rgi_path = "" # Initialize RGI path to be accessible outside the try block
    try
        rgi_path = joinpath(prepro_dir, simparams.rgi_paths[rgi_id])
    catch
        @error "RGI path not found for: $rgi_id"
    end

    period = _climate_period(simparams)
    raw_climate_clipped_file = "raw_climate_$(period).nc"
    if !isfile(joinpath(rgi_path, raw_climate_clipped_file))
        println("Getting raw climate data for: ", rgi_id)
        # Get raw climate data for gdir. We start a year before the period to ensure we have enough data
        # for variables with a sliding time window
        tspan_date = partial_year(Day, period[1] - 1):Day(1):partial_year(
            Day, period[2])
        climate = get_raw_climate_data(rgi_path, simparams.climate_data_source)
        # Make sure the desired period is covered by the climate data
        period = trim_period(tspan_date, climate)
        climTstart = Date(dims(climate, Ti)[begin])
        climTend = Date(dims(climate, Ti)[end])

        if any((climTstart <= period[begin]) & any(climTend >= period[end]))
            climate = _slice_climate_between_dates(climate, period[begin], period[end])
        else
            slice_start = max(period[begin], climTstart)
            slice_end = min(period[end], climTend)
            if slice_start <= slice_end
                @warn "No overlapping period available between climate tspan! Climate data range from $(climTstart) to $(climTend)."
                climate = _slice_climate_between_dates(climate, slice_start, slice_end)
            else
                error("No climate data available for the requested period $(period[begin]) to $(period[end]). Climate data range from $(climTstart) to $(climTend).")
            end
        end
        # Save raw gdir climate on disk
        write(joinpath(rgi_path, raw_climate_clipped_file), climate)
        GC.gc()
    end
end

"""
    get_cumulative_climate!(climate, t::Sleipnir.Float, step::Sleipnir.Float, gradient_bounds=[-0.009, -0.003])
    get_cumulative_climate!(climate, period::StepRange{Date, Day}, gradient_bounds=[-0.009, -0.003])

Calculate and update the cumulative climate data for a given period.
The user can choose between providing a specific time `t` and a time step `step`, or a time period defined by `period`.

# Keyword arguments

  - `climate::Climate`: The climate object containing raw climate data.
  - `gradient_bounds::Vector{Float64}`: Optional. The bounds within which to clamp the gradient values. Default is `[-0.009, -0.003]`.
    Optional parameters to specify the time period:
  - `t::Sleipnir.Float`: Time at which the cumulative climate data should be computed.
  - `step::Sleipnir.Float`: Time step used to compute the cumulative climate data. Together with `t` they define a time period.
    or
  - `period::StepRange{Date, Day}`: The time period for which to compute the cumulative climate data.

# Updates

  - `climate.climate_raw_step`: The raw climate data for the given period.
  - `climate.avg_temps`: The average temperature for the given period.
  - `climate.avg_gradients`: The average gradient for the given period.
  - `climate.climate_step.prcp`: The cumulative precipitation for the given period.
  - `climate.climate_step.temp`: The cumulative temperature for the given period.
  - `climate.climate_step.gradient`: The cumulative gradient for the given period.
  - `climate.climate_step.avg_temp`: The average temperature for the given period.
  - `climate.climate_step.avg_gradient`: The average gradient for the given period.
  - `climate.climate_step.ref_hgt`: The reference height from the raw climate data.
"""
function get_cumulative_climate!(
        climate, t::Sleipnir.Float, step::Sleipnir.Float, gradient_bounds = [
            -0.009, -0.003])
    # First we get the dates of the current time and the previous step
    period = partial_year(Day, t - step):Day(1):partial_year(Day, t)
    get_cumulative_climate!(climate, period, gradient_bounds)
end

function get_cumulative_climate!(
        climate, period::StepRange{Date, Day}, gradient_bounds = [-0.009, -0.003])
    climate.climate_raw_step = _slice_climate_between_dates(
        climate.raw_climate, period[begin], period[end])

    climate.avg_temps = mean(climate.climate_raw_step.temp)

    climate.avg_gradients = mean(climate.climate_raw_step.gradient)
    climate.climate_raw_step.gradient.data .= clamp.(
        climate.climate_raw_step.gradient.data, gradient_bounds[1], gradient_bounds[2]) # Clip gradients within plausible values
    climate.climate_step.prcp = round(sum(climate.climate_raw_step.prcp); digits = 8)
    # Sum daily PDDs without overwriting temp so original per-day temperatures
    # remain available for the daily snow/PDD downscaling in downscale_2D_climate!.
    climate.climate_step.temp = round(
        sum(max.(climate.climate_raw_step.temp.data, 0)); digits = 8)
    climate.climate_step.gradient = round(
        sum(climate.climate_raw_step.gradient); digits = 8)
    climate.climate_step.albedo = _aggregate_raw_layer(
        climate.climate_raw_step, :fal; reducer = mean)
    climate.climate_step.slhf = _aggregate_raw_layer(climate.climate_raw_step, :slhf)
    climate.climate_step.sshf = _aggregate_raw_layer(climate.climate_raw_step, :sshf)
    climate.climate_step.ssrd = _aggregate_raw_layer(climate.climate_raw_step, :ssrd)
    climate.climate_step.str = _aggregate_raw_layer(climate.climate_raw_step, :str)
    climate.climate_step.avg_temp = round(climate.avg_temps; digits = 8)
    climate.climate_step.avg_gradient = round(climate.avg_gradients; digits = 8)
    climate.climate_step.ref_hgt = round(climate.ref_hgt; digits = 8)
end

"""
    get_cumulative_climate(
        climate::RasterStack,
        gradient_bounds::Vector{Float64}=[-0.009, -0.003],
    )

Calculate cumulative climate statistics from the given climate data.

# Keyword arguments

  - `climate::RasterStack`: A `RasterStack` object containing temperature, precipitation, and gradient data.
  - `gradient_bounds::Vector{Float64}`: A two-element vector specifying the lower and upper bounds for the gradient values. Defaults to `[-0.009, -0.003]`.

# Returns

  - `climate_sum::ClimateStep`: A struct containing the following fields:

      + `"temp"`: The sum of positive degree days (PDDs) from the temperature data.
      + `"prcp"`: The sum of precipitation data.
      + `"gradient"`: The sum of gradient data, clipped within the specified bounds.
      + `"avg_temp"`: The average temperature.
      + `"avg_gradient"`: The average gradient.
      + `"ref_hgt"`: The reference height from the climate metadata.

# Notes

  - The temperature data is modified to only include positive degree-day values (PDDs).
  - The gradient data is clipped within the specified bounds to ensure plausible values.
"""
function get_cumulative_climate(
        climate::RasterStack,
        gradient_bounds::Vector{Float64} = [-0.009, -0.003]
)
    avg_temp = mean(climate.temp)
    avg_gradient = mean(climate.gradient)
    copy_climate = deepcopy(climate)
    copy_climate.temp.data .= max.(copy_climate.temp.data, 0.0) # get PDDs
    copy_climate.gradient.data .= clamp.(
        copy_climate.gradient.data, gradient_bounds[1], gradient_bounds[2]) # Clip gradients within plausible values
    climate_sum = ClimateStep(
        temp = round(sum(copy_climate.temp); digits = 8),
        prcp = round(sum(climate.prcp); digits = 8),
        gradient = round(sum(copy_climate.gradient); digits = 8),
        albedo = _aggregate_raw_layer(climate, :fal; reducer = mean),
        slhf = _aggregate_raw_layer(climate, :slhf),
        sshf = _aggregate_raw_layer(climate, :sshf),
        ssrd = _aggregate_raw_layer(climate, :ssrd),
        str = _aggregate_raw_layer(climate, :str),
        avg_temp = round(avg_temp; digits = 8),
        avg_gradient = round(avg_gradient; digits = 8),
        ref_hgt = round(Sleipnir.Float(metadata(climate)["ref_hgt"]); digits = 8)
    )
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
function get_raw_climate_data(rgi_path::String, climate_data_source::Symbol)
    if climate_data_source == :W5E5
        climate = RasterStack(joinpath(rgi_path, "climate_historical_daily_W5E5.nc"))
    elseif climate_data_source == :ERA5
        monthly_path = joinpath(rgi_path, "climate_historical_monthly_ERA5.nc")
        daily_path = joinpath(rgi_path, "climate_historical_daily_ERA5.nc")
        if isfile(monthly_path)
            climate = RasterStack(monthly_path)
        elseif isfile(daily_path)
            climate = RasterStack(daily_path)
        else
            throw(ArgumentError(
                "No ERA5 climate file found in $(rgi_path). Expected climate_historical_monthly_ERA5.nc or climate_historical_daily_ERA5.nc."
            ))
        end
    else
        throw(ArgumentError("Unsupported climate data source"))
    end
    return climate
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
    climate.temp.data .= climate.temp.data .+
                         climate.gradient.data .*
                         (mean(dem.data[:]) .- metadata(climate)["ref_hgt"])
end

function apply_t_grad!(climate::RasterStack, dem::Matrix{<: AbstractFloat})
    # We apply the gradients to the temperature
    climate.temp.data .= climate.temp.data .+
                         climate.gradient.data .*
                         (mean(dem) .- metadata(climate)["ref_hgt"])
end

function apply_t_grad_gridded(climate::RasterStack, dem::Matrix{<: AbstractFloat})
    # We apply the gradients to the temperature
    dummy_grid = zeros(size(dem))
    temps_2D = [temp .+ dummy_grid for temp in climate.temp.data]
    for i in eachindex(temps_2D)
        temps_2D[i] .= temps_2D[i] .+
                       climate.gradient.data[i] .* (dem .- metadata(climate)["ref_hgt"])
    end

    return temps_2D
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

    # Update 2D climate structure

  - The function modifies the `glacier` object in place.
"""
function downscale_2D_climate!(
        glacier::Glacier2D;
        include_topography::Bool = false,
        topography_window_m::Sleipnir.Float = Sleipnir.Float(200.0),
        temp_bias::Sleipnir.Float = Sleipnir.Float(0.0))
    climate = glacier.climate
    @. climate.climate_2D_step.elevation_diff = glacier.S - climate.climate_step.ref_hgt
    if include_topography
        climate.climate_2D_step.slope,
        climate.climate_2D_step.aspect = compute_surface_topography(
            glacier; window_m = topography_window_m)
    end
    climate.climate_2D_step.albedo .= climate.climate_step.albedo
    climate.climate_2D_step.slhf .= climate.climate_step.slhf
    climate.climate_2D_step.sshf .= climate.climate_step.sshf
    climate.climate_2D_step.ssrd .= climate.climate_step.ssrd
    climate.climate_2D_step.str .= climate.climate_step.str
    climate.climate_2D_step.gradient = climate.climate_step.gradient
    climate.climate_2D_step.avg_gradient = climate.climate_step.avg_gradient

    # Per-day snow/PDD with per-day elevation gradient (daily TI model).
    # climate_raw_step.temp holds original daily temperatures; gradient is clamped.
    temp_vec = vec(climate.climate_raw_step.temp.data)
    prcp_vec = vec(climate.climate_raw_step.prcp.data)
    grad_vec = vec(climate.climate_raw_step.gradient.data)
    ΔS = climate.climate_2D_step.elevation_diff
    fill!(climate.climate_2D_step.snow, 0)
    fill!(climate.climate_2D_step.PDD, 0)
    for d in eachindex(temp_vec)
        T_d, g_d, p_d = temp_vec[d], grad_vec[d], prcp_vec[d]
        T_eff = T_d + temp_bias
        @. climate.climate_2D_step.snow += p_d * clamp((2 - T_eff - g_d * ΔS) / 2, 0, 1)
        @. climate.climate_2D_step.PDD += max(0, T_eff + g_d * ΔS)
    end
    climate.climate_2D_step.rain .= climate.climate_step.prcp .-
                                    climate.climate_2D_step.snow
    climate.climate_2D_step.temp .= climate.climate_step.avg_temp .+ temp_bias .+
                                    climate.climate_2D_step.avg_gradient .* ΔS
end

"""
    downscale_2D_climate(climate_step::ClimateStep, S::Matrix{<: AbstractFloat}, Coords::Dict)

Downscales climate data to a 2D grid based on the provided matrix of surface elevation and coordinates.

# Arguments

  - `climate_step::ClimateStep`: A struct containing climate data for a specific time step. Expected fields are:

      + `"avg_temp"`: Average temperature.
      + `"temp"`: Temperature.
      + `"prcp"`: Precipitation.
      + `"gradient"`: Temperature gradient.
      + `"avg_gradient"`: Average temperature gradient.
      + `"ref_hgt"`: Reference height.

  - `S::Matrix{<: AbstractFloat}`: Surface elevation data.

  - `Coords::Dict`: A dictionary with keys `"lon"` and `"lat"` for longitude and latitude coordinates.

# Returns

  - `Climate2Dstep`: A `Climate2Dstep` object containing the downscaled climate data with fields:

      + `temp`: 2D array of temperature.
      + `PDD`: 2D array of positive degree days.
      + `snow`: 2D array of snow precipitation.
      + `rain`: 2D array of rain precipitation.
      + `gradient`: Temperature gradient.
      + `avg_gradient`: Average temperature gradient.
      + `x`: Longitude coordinates.
      + `y`: Latitude coordinates.
      + `ref_hgt`: Reference height.

# Description    # Create dummy 2D arrays to have a base to apply gradients afterwards

This function creates dummy 2D arrays based on the provided surface elevation data and applies the climate step data to these arrays. It then constructs a `Climate2Dstep` object with the downscaled climate data and applies temperature gradients to compute the snow/rain fraction for the selected period.
"""
function downscale_2D_climate(
        climate_step::ClimateStep,
        climate_raw_step,
        S::Matrix{<: AbstractFloat},
        Coords::Dict;
        include_topography::Bool = false,
        topography_window_m::Sleipnir.Float = Sleipnir.Float(200.0),
        Δx::Union{Nothing, Sleipnir.Float} = nothing,
        Δy::Union{Nothing, Sleipnir.Float} = nothing,
        temp_bias::Sleipnir.Float = Sleipnir.Float(0.0))
    dummy_grid = zeros(size(S))
    temp_2D = climate_step.avg_temp .+ dummy_grid
    elevation_diff_2D = S .- climate_step.ref_hgt
    PDD_2D = zeros(Sleipnir.Float, size(S))
    snow_2D = zeros(Sleipnir.Float, size(S))
    rain_2D = zeros(Sleipnir.Float, size(S))
    slope_2D = zero(Sleipnir.Float) .+ dummy_grid
    aspect_2D = zero(Sleipnir.Float) .+ dummy_grid
    if include_topography
        Δx === nothing && error("Missing Δx for topography computation")
        Δy === nothing && error("Missing Δy for topography computation")
        slope_2D,
        aspect_2D = compute_surface_topography(
            S,
            Δx,
            Δy;
            window_m = topography_window_m)
    end
    albedo_2D = zero(Sleipnir.Float) .+ dummy_grid
    slhf_2D = zero(Sleipnir.Float) .+ dummy_grid
    sshf_2D = zero(Sleipnir.Float) .+ dummy_grid
    ssrd_2D = zero(Sleipnir.Float) .+ dummy_grid
    str_2D = zero(Sleipnir.Float) .+ dummy_grid

    albedo_2D .= climate_step.albedo
    slhf_2D .= climate_step.slhf
    sshf_2D .= climate_step.sshf
    ssrd_2D .= climate_step.ssrd
    str_2D .= climate_step.str

    climate_2D_step = Climate2Dstep{Sleipnir.Float}(
        temp = temp_2D,
        PDD = PDD_2D,
        snow = snow_2D,
        rain = rain_2D,
        elevation_diff = elevation_diff_2D,
        aspect = aspect_2D,
        albedo = albedo_2D,
        slhf = slhf_2D,
        slope = slope_2D,
        sshf = sshf_2D,
        ssrd = ssrd_2D,
        str = str_2D,
        gradient = Float64(climate_step.gradient),
        avg_gradient = Float64(climate_step.avg_gradient),
        x = Coords["lon"],
        y = Coords["lat"],
        ref_hgt = Float64(climate_step.ref_hgt)
    )

    temp_vec = vec(climate_raw_step.temp.data)
    prcp_vec = vec(climate_raw_step.prcp.data)
    grad_vec = vec(climate_raw_step.gradient.data)
    ΔS = elevation_diff_2D
    for d in eachindex(temp_vec)
        T_d, g_d, p_d = temp_vec[d], grad_vec[d], prcp_vec[d]
        T_eff = T_d + temp_bias
        @. climate_2D_step.snow += p_d * clamp((2 - T_eff - g_d * ΔS) / 2, 0, 1)
        @. climate_2D_step.PDD += max(0, T_eff + g_d * ΔS)
    end
    climate_2D_step.rain .= climate_step.prcp .- climate_2D_step.snow
    climate_2D_step.temp .+= temp_bias .+ climate_2D_step.avg_gradient .* ΔS

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
  - `float::Sleipnir.Float`: The floating-point year value.

# Returns

  - `Date`: The calculated date corresponding to the partial year.
"""
function partial_year(period::Type{<:Period}, float::Sleipnir.Float)
    _year, Δ = divrem(float, 1)
    year_start = Date(convert(Int, _year))
    year = period((year_start + Year(1)) - year_start)
    partial = period(round(Dates.value(year) * Δ))
    year_start + partial
end
function partial_year(period::Type{<:Period}, floats::Vector{Sleipnir.Float})
    map(f -> partial_year(period, f), floats)
end

"""
    partial_year(float::Sleipnir.Float) -> Sleipnir.Float

Calculate the partial year value based on the given floating-point number.

# Arguments

  - `float::Sleipnir.Float`: A floating-point number representing the fraction of the year.

# Returns

  - `Sleipnir.Float`: The calculated partial year value.
"""
partial_year(float) = partial_year(Day, float)

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
function get_longterm_temps(rgi_id::String, params::Parameters,
        climate::RasterStack, S::Matrix{<: AbstractFloat})
    glacier_gd = RasterStack(joinpath(
        prepro_dir, params.simulation.rgi_paths[rgi_id], "gridded_data.nc"))
    temp_orig = copy(climate.temp.data)  # avoid corrupting temp in-place
    apply_t_grad!(climate, glacier_gd.topo)
    temps_2D = apply_t_grad_gridded(climate, S)

    # Scalar long-term temps
    longterm_temps_scalar = mean.(groupby(climate.temp, Ti=>year)).data
    # Gridded long-term temps
    longterm_temps_gridded = mean(temps_2D, dims = 1)[1]

    climate.temp.data .= temp_orig
    return longterm_temps_scalar, longterm_temps_gridded
end

"""
    get_winter_prcp_factor(glacier::AbstractGlacier, params::Parameters; prcp_fac_bounds) -> Float

Compute a glacier-specific precipitation correction factor from mean winter precipitation,
following OGGM's `decide_winter_precip_factor`.

# Arguments

  - `glacier::AbstractGlacier`: Glacier providing `rgi_id`, `cenlat` (hemisphere) and
    `climate.climate_data_source` (`:W5E5` or `:ERA5`).
  - `params::Parameters`: Simulation parameters, used to locate the climate file.
  - `prcp_fac_bounds`: Lower/upper clip bounds (default `(0.1, 10.0)`, matching OGGM).

# Description

Reads the full daily climate record (1979–2019, the window OGGM's coefficients were fit on),
averages winter daily precipitation (NH: Oct–Apr, SH: Apr–Oct) in kg/m²/day, and maps it to a
precipitation factor: log fit for W5E5, linear fit for ERA5. If the daily file is missing, falls
back to 2.5. Note: we average daily values directly (OGGM averages monthly daily-rates; equivalent
within noise).
"""
function get_winter_prcp_factor(glacier::AbstractGlacier, params::Parameters;
        prcp_fac_bounds::Tuple{<: Real, <: Real} = (0.1, 10.0))
    source = glacier.climate.climate_data_source
    fpath = joinpath(prepro_dir, params.simulation.rgi_paths[glacier.rgi_id],
        "climate_historical_daily_$(source).nc")
    if !isfile(fpath)
        @warn "get_winter_prcp_factor: daily climate file not found for $(glacier.rgi_id) " *
              "($(fpath)). Falling back to prcp_fac = 2.5."
        return Sleipnir.Float(2.5)
    end
    clim = RasterStack(fpath)

    # Winter months by hemisphere; restrict to the 1979–2019 coefficient-fit window
    winter = glacier.cenlat >= 0 ? (10, 11, 12, 1, 2, 3, 4) : (4, 5, 6, 7, 8, 9, 10)
    dates = Date.(collect(dims(clim, Ti)))
    sel = [(month(d) in winter) && (Date(1979, 1, 1) <= d <= Date(2019, 12, 31))
           for d in dates]
    w_prcp = mean(clim.prcp.data[sel])  # kg/m²/day

    # OGGM winter-precip relationships (Nov 2025 refit; t_melt=-1, cte lapse, monthly)
    prcp_fac = if source == :W5E5
        -1.0614 * log(w_prcp) + 3.9200          # log fit
    elseif source == :ERA5
        -0.09078476 * w_prcp + 2.43505368       # linear fit
    else
        throw(ArgumentError(
            "get_winter_prcp_factor: unsupported climate_data_source $(source)."))
    end

    return Sleipnir.Float(clamp(prcp_fac, prcp_fac_bounds[1], prcp_fac_bounds[2]))
end
