
###############################################
############  FUNCTIONS   #####################
###############################################

export initialize_glacier_climate!, downscale_2D_climate!, downscale_2D_climate,
        get_cumulative_climate!, get_cumulative_climate, apply_t_cumul_grad!,
         apply_t_grad!, trim_period, partial_year, get_longterm_temps

using Dates # to provide correct Julian time slices

"""
    function initialize_glacier_climate!(glacier::AbstractGlacier, params::Parameters)

Initializes the `Climate` data structure for a given `Glacier``
"""
function initialize_glacier_climate!(glacier::AbstractGlacier, params::Parameters)

    dummy_period = partial_year(Day, params.simulation.tspan[1]):Day(1):partial_year(Day, params.simulation.tspan[1] + params.simulation.step)
    raw_climate = RasterStack(joinpath(params.simulation.rgi_paths[glacier.rgi_id], "raw_climate_$(params.simulation.tspan).nc"))
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

function generate_raw_climate_files(rgi_id::String, simparams::SimulationParameters) where {F <: AbstractFloat}
    rgi_path = simparams.rgi_paths[rgi_id]
    if !ispath(joinpath(rgi_path, "raw_climate_$(simparams.tspan).nc"))
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
    get_cumulative_climate(climate, gradient_bounds=[-0.009, -0.003], default_grad=-0.0065)

Computes Positive Degree Days (PDDs) and cumulative rainfall and snowfall from climate data.
"""
function get_cumulative_climate!(climate, period, gradient_bounds=[-0.009, -0.003], default_grad=-0.0065)
    climate.climate_raw_step[] = climate.raw_climate.sel(time=period)
    climate.avg_temps[] = climate.climate_raw_step[].temp.mean() 

    climate.avg_gradients[] = climate.climate_raw_step[].gradient.mean() 
    climate.climate_raw_step[].temp.data = climate.climate_raw_step[].temp.where(climate.climate_raw_step[].temp > 0.0, 0.0).data # get PDDs
    climate.climate_raw_step[].gradient.data = utils[].clip_array(climate.climate_raw_step[].gradient.data, gradient_bounds[1], gradient_bounds[2]) # Clip gradients within plausible values
    climate.climate_step[] = climate.climate_raw_step[].sum() # get monthly cumulative values
    climate.climate_step[] = climate.climate_step[].assign(Dict("avg_temp"=>climate.avg_temps[])) 
    climate.climate_step[] = climate.climate_step[].assign(Dict("avg_gradient"=>climate.avg_gradients[]))
    climate.climate_step[].attrs = climate.climate_raw_step[].attrs
end

function get_cumulative_climate(climate, gradient_bounds=[-0.009, -0.003], default_grad=-0.0065)
    climate.temp.data .= max.(climate.temp.data, 0.0) # get PDDs
    climate.gradient.data .= clamp.(climate.gradient.data, gradient_bounds[1], gradient_bounds[2]) # Clip gradients within plausible values attributes 
    climate_sum = Dict("temp" => sum(climate.temp), 
                       "prcp" => sum(climate.prcp), 
                       "avg_temp" => mean(climate.temp), 
                       "avg_gradient" => mean(climate.gradient),
                       "ref_hgt" => metadata(climate)["ref_hgt"])
    return climate_sum
end

"""
    get_raw_climate_data(rgi_path::String)

Load the netCDF file containing the climate data for that glacier.
"""
function get_raw_climate_data(rgi_path::String)
    climate = RasterStack(joinpath(rgi_path, "climate_historical_daily_W5E5.nc"))
    return climate
end

# TODO: make snow/rain thresholds customizable 
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
    apply_t_grad!(climate, g_dem)

Applies temperature gradients to the glacier 2D climate data based on a DEM.  
"""
function apply_t_grad!(climate::RasterStack, dem::Raster)
    # We apply the gradients to the temperature
    climate.temp.data .= climate.temp.data .+ climate.gradient.data .* (mean(dem.data[:]) .- climate.ref_hgt) 
end

"""
    downscale_2D_climate(climate, g_dem)

Projects climate data to the glacier matrix by simply copying the closest gridpoint to all matrix gridpoints.
Generates a new xarray Dataset which is returned.   
"""
function downscale_2D_climate!(glacier::Glacier2D)
    # Update 2D climate structure
    climate = glacier.climate
    climate.climate_2D_step.temp .= climate.climate_step.avg_temp.data
    climate.climate_2D_step.PDD .= climate.climate_step.temp.data
    climate.climate_2D_step.snow .= climate.climate_step.prcp.data
    climate.climate_2D_step.rain .= climate.climate_step.prcp.data
    # Update gradients
    climate.climate_2D_step.gradient .= climate.climate_step.gradient.data
    climate.climate_2D_step.avg_gradient .= climate.climate_step.avg_gradient.data

    # Apply temperature gradients and compute snow/rain fraction for the selected period
    apply_t_cumul_grad!(climate.climate_2D_step, reshape(glacier.S, size(glacier.S))) # Reproject current S with xarray structure
end

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
                       gradient=climate_step["gradient"],
                       avg_gradient=climate_step["avg_gradient"],
                       x=glacier.S_coords.x.data,
                       y=glacier.S_coords.y.data,
                       ref_hgt=climate_step["ref_hgt"])

    # Apply temperature gradients and compute snow/rain fraction for the selected period
    apply_t_cumul_grad!(climate_2D_step, reshape(glacier.S, size(glacier.S))) # Reproject current S with xarray structure

    return climate_2D_step

end

function downscale_2D_climate(glacier::Glacier2D)
    climate_2D_step = downscale_2D_climate(glacier.climate.climate_step[], glacier)
    return climate_2D_step
end

"""
    trim_period(period, climate)

Trims a time period based on the time range of a climate series. 
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

function partial_year(period::Type{<:Period}, float)
    _year, Δ = divrem(float, 1)
    year_start = Date(_year)
    year = period((year_start + Year(1)) - year_start)
    partial = period(round(Dates.value(year) * Δ))
    year_start + partial
end
partial_year(float) = partial_year(Day, float)


function get_longterm_temps(rgi_id::String, params::Parameters)
    rgi_path = params.simulation.rgi_paths[rgi_id]
    glacier_gd = RasterStack(joinpath(rgi_path, "gridded_data.nc"))
    climate = RasterStack(joinpath(rgi_path, "raw_climate_$(params.simulation.tspan).nc"))
    apply_t_grad!(climate, glacier_gd.topo)
    longterm_temps = mean.(groupby(climate.temp, Ti=>year)).data
    return longterm_temps
end

function get_longterm_temps(rgi_id::String, params::Parameters, climate::RasterStack)
    glacier_gd = RasterStack(joinpath(params.simulation.rgi_paths[rgi_id], "gridded_data.nc"))
    apply_t_grad!(climate, glacier_gd.topo)
    longterm_temps = mean.(groupby(climate.temp, Ti=>year)).data
    return longterm_temps
end


