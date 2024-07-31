
###############################################
############  FUNCTIONS   #####################
###############################################

export initialize_glacier_climate!, downscale_2D_climate!, downscale_2D_climate,
        get_cumulative_climate!, get_cumulative_climate, apply_t_cumul_grad!,
         apply_t_grad!, trim_period, partial_year, get_longterm_temps

using Dates # to provide correct Julian time slices 

"""
    function initialize_glacier_climate!(glacier::Glacier, params::Parameters)

Initializes the `Climate` data structure for a given `Glacier``
"""
function initialize_glacier_climate!(glacier::AbstractGlacier, params::Parameters)

    dummy_period = partial_year(Day, params.simulation.tspan[1]):Day(1):partial_year(Day, params.simulation.tspan[1] + params.simulation.step)
    raw_climate = xr[].open_dataset(joinpath(pyconvert(String,glacier.gdir.dir), "raw_climate_$(params.simulation.tspan).nc"))
    climate_step = Ref{Py}(get_cumulative_climate(raw_climate.sel(time=dummy_period)))
    climate_2D_step = downscale_2D_climate(climate_step[], glacier)
    longterm_temps = get_longterm_temps(glacier.gdir, raw_climate)
    glacier.climate = Climate2D(raw_climate = raw_climate,
                            climate_raw_step = Ref{Py}(raw_climate.sel(time=dummy_period)),
                            #climate_cum_step = raw_climate.sel(time=dummy_period).sum(),
                            climate_step = climate_step,
                            climate_2D_step = climate_2D_step,
                            longterm_temps = pyconvert(Vector,longterm_temps),
                            avg_temps = Ref{Py}(raw_climate.sel(time=dummy_period).temp.mean()),
                            avg_gradients = Ref{Py}(raw_climate.sel(time=dummy_period).gradient.mean()))

end

function generate_raw_climate_files(gdir::Py, tspan::Tuple{F, F}) where {F <: AbstractFloat}
    if !ispath(joinpath(pyconvert(String,gdir.dir), "raw_climate_$tspan.nc"))
        println("Getting raw climate data for: ", gdir.rgi_id)
        # Get raw climate data for gdir
        tspan_date = partial_year(Day, tspan[1]):Day(1):partial_year(Day, tspan[2])
        climate =  get_raw_climate_data(gdir)
        # Make sure the desired period is covered by the climate data
        period = trim_period(tspan_date, climate) 
        if any((jldate(climate.time, 0) <= period[1]) & any(jldate(climate.time, -1) >= period[end]))
            climate = climate.sel(time=period) # Crop desired time period
        else
            @warn "No overlapping period available between climate tspan!" 
        end
        # Save raw gdir climate on disk 
        climate.to_netcdf(joinpath(pyconvert(String,gdir.dir), "raw_climate_$tspan.nc"))
        climate.close()
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
    avg_temp = climate.temp.mean() 
    avg_gradients = climate.gradient.mean() 
    climate.temp.data = climate.temp.where(climate.temp > 0, 0).data # get PDDs
    climate.gradient.data = utils[].clip_array(climate.gradient.data, gradient_bounds[1], gradient_bounds[2]) # Clip gradients within plausible values
    attributes = climate.attrs
    climate_sum = climate.sum() # get monthly cumulative values
    climate_sum = climate_sum.assign(Dict("avg_temp"=>avg_temp)) 
    climate_sum = climate_sum.assign(Dict("avg_gradient"=>avg_gradients))
    climate_sum.attrs = attributes
    return climate_sum
end

"""
    get_raw_climate_data(gdir, temp_resolution="daily", climate="W5E5")

Downloads the raw W5E5 climate data with a given resolution (daily by default). Returns an xarray Dataset. 
"""
function get_raw_climate_data(gdir; temp_resolution="daily", climate="W5E5")
    MBsandbox[].process_w5e5_data(gdir, climate_type=climate, temporal_resol=temp_resolution) 
    fpath = gdir.get_filepath("climate_historical", filesuffix="_daily_W5E5")
    climate = xr[].open_dataset(fpath)
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
function apply_t_grad!(climate::Py, dem)
    # We apply the gradients to the temperature
    # /!\ AVOID USING `.=` IN JULIA TO ASSIGN. IT'S NOT HANDLED BY XARRAY. USE `=` INSTEAD
    climate.temp.data = climate.temp.data + climate.gradient.data .* (mean(dem.data.flatten()) .- climate.ref_hgt)
end

"""
    downscale_2D_climate(climate, g_dem)

Projects climate data to the glacier matrix by simply copying the closest gridpoint to all matrix gridpoints.
Generates a new xarray Dataset which is returned.   
"""
function downscale_2D_climate!(glacier::Glacier2D)
    # Update 2D climate structure
    climate = glacier.climate
    FT = eltype(glacier.S[1])

    climate.climate_2D_step.temp .= pyconvert(FT,climate.climate_step[].avg_temp.data[()])
    climate.climate_2D_step.PDD .= pyconvert(FT,climate.climate_step[].temp.data[()])
    climate.climate_2D_step.snow .= pyconvert(FT,climate.climate_step[].prcp.data[()])
    climate.climate_2D_step.rain .= pyconvert(FT,climate.climate_step[].prcp.data[()]) 
    # Update gradients
    climate.climate_2D_step.gradient[] = pyconvert(FT,climate.climate_step[].gradient.data[()])
    climate.climate_2D_step.avg_gradient[] = pyconvert(FT,climate.climate_step[].avg_gradient.data[()])

    # Apply temperature gradients and compute snow/rain fraction for the selected period
    apply_t_cumul_grad!(climate.climate_2D_step, reshape(glacier.S, size(glacier.S))) # Reproject current S with xarray structure
end

function downscale_2D_climate(climate_step::Py, glacier::Glacier2D)
    # Create dummy 2D arrays to have a base to apply gradients afterwards
    FT = typeof(glacier.S[1])
    dummy_grid = ones(size(glacier.S))
    temp_2D = pyconvert(FT, climate_step.avg_temp.data[()]) .* dummy_grid
    PDD_2D = pyconvert(FT, climate_step.temp.data[()]) .* dummy_grid
    snow_2D = pyconvert(FT, climate_step.prcp.data[()]) .* dummy_grid
    rain_2D = pyconvert(FT, climate_step.prcp.data[()]) .* dummy_grid

    # We generate a new dataset with the scaled data
    climate_2D_step = Climate2Dstep(temp=temp_2D,
                                   PDD=PDD_2D,
                                   snow=snow_2D,
                                   rain=rain_2D,
                                   gradient=Ref{FT}(pyconvert(FT,climate_step.gradient.data[()])),
                                   avg_gradient=Ref{FT}(pyconvert(FT,climate_step.avg_gradient.data[()])),
                                   x=pyconvert(Vector{FT},glacier.S_coords.x.data),
                                   y=pyconvert(Vector{FT},glacier.S_coords.y.data),
                                   ref_hgt=Ref{FT}(pyconvert(FT,climate_step.ref_hgt)))

    # Apply temperature gradients and compute snow/rain fraction for the selected period
    apply_t_cumul_grad!(climate_2D_step, reshape(glacier.S, size(glacier.S))) # Reproject current S with xarray structure
    
    return climate_2D_step

end

function downscale_2D_climate(glacier::Glacier2D)
    climate_2D_step = downscale_2D_climate(glacier.climate.climate_step[], glacier)
    return climate_2D_step
end

"""
    jldate(pydate, 0)

Converts a Python date (generally from xarray) to a Julia `Date`. WARNING: it requires Python indices, e.g. 
0 for the beginning and -1 for the end.  
"""
function jldate(pydate, idx)
    return Date(pyconvert(Int, pydate.dt.year.data[idx]), pyconvert(Int,pydate.dt.month.data[idx]), pyconvert(Int,pydate.dt.day.data[idx]))
end

"""
    trim_period(period, climate)

Trims a time period based on the time range of a climate series. 
"""
function trim_period(period, climate)
    if any(jldate(climate.time, 0) > period[1])
        head = jldate(climate.time, 0)
        period = Date(year(head), 10, 1):Day(1):period[end] # make it a hydrological year
    end
    if any(jldate(climate.time, -1) > period[end])
        tail = jldate(climate.time, -1)
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


function get_longterm_temps(gdir::Py, tspan)
    climate = xr[].open_dataset(joinpath(gdir.dir, "raw_climate_$tspan.nc")) # load only once at the beginning
    dem = rioxarray[].open_rasterio(gdir.get_filepath("dem"))
    apply_t_grad!(climate, dem)
    longterm_temps = climate.groupby("time.year").mean().temp.data
    return longterm_temps
end

function get_longterm_temps(gdir::Py, climate::Py)
    dem = rioxarray[].open_rasterio(gdir.get_filepath("dem"))
    apply_t_grad!(climate, dem)
    longterm_temps = climate.groupby("time.year").mean().temp.data
    return longterm_temps
end


