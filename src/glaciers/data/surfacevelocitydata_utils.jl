export initialize_surfacevelocitydata
export plot_timeseries, plot_count

"""
    initialize_surfacevelocitydata(file::String; interp=false)

Initialize SurfaceVelocityData from Rabatel et. al (2023). 

Arguments
=================
    - `file`: name of netcdf file with data
    - `interp`: boolean variable indicating if we use the inporpolated data or not
"""
function initialize_surfacevelocitydata(file::String; interp=false)

    # Date of first adquisition
    date1 = ncread(file, "date1")
    # Date of second adquisition
    date2 = ncread(file, "date2")

    # Middle date
    if !interp 
        date_mean = 0.5 .* date1 .+ 0.5 .* date2
        date_mean_offset = date1_offset = date2_offset = 0
    else 
        date_mean = ncread(file, "mid_date")
        # Extract string from metadata of when days start counting
        # Default format of attribute "units" is "days since 2015-08-14 00:00:00"
        date_mean_since = ncgetatt(file, "mid_date", "units")[12:21] # e.g., "2015-08-14"
        date1_offset_since = ncgetatt(file, "date1", "units")[12:21] # e.g., "2015-07-30"
        date2_offset_since = ncgetatt(file, "date2", "units")[12:21] # e.g., "2015-08-29"
        # Convertion to Julia datetime
        date_mean_offset = datetime2julian(DateTime(date_mean_since)) - 2400000.5
        date1_offset = datetime2julian(DateTime(date1_offset_since)) - 2400000.5
        date2_offset = datetime2julian(DateTime(date_mean_offset)) - 2400000.5
    end

    # Convert dates from Modified Julian Days to datetipes
    date1 = mjd.(date1 .+ date1_offset)
    date2 = mjd.(date2 .+ date2_offset)
    date_mean = mjd.(date_mean .+ date_mean_offset)
    date_error = date2 .- date1

    # We convert from Date to DateTime
    date_mean = DateTime.(date_mean)
    date1 = DateTime.(date1)
    date2 = DateTime.(date2)

    # Read data from netcdf file
    x = ncread(file, "x")
    y = ncread(file, "y")

    # Velocity in the x direction (m/yr)
    vx = ncread(file, "vx")              
    vy = ncread(file, "vy")    

    # Compute absolute velocity
    vabs = (vx.^2 .+ vy.^2).^0.5  
    # The sqrt operation in Julia promotes Float32 to Float64. We convert manually 
    # to keep type consistency     
    vabs = convert(typeof(vx), vabs)   
    
    # Error is reported once per timespan, so upper bounds are given by absolute error
    if !interp
        vx_error = ncread(file, "error_vx")
        vy_error = ncread(file, "error_vy")
        # Absolute error uncertanty using propagation of uncertanties 
        vx_ratio_max = map(i -> max_or_empty(abs.(vx[:,:,i][vabs[:,:,i] .> 0.0]) ./ vabs[:,:,i][vabs[:,:,i] .> 0.0]), 1:size(vx)[3])
        vy_ratio_max = map(i -> max_or_empty(abs.(vy[:,:,i][vabs[:,:,i] .> 0.0]) ./ vabs[:,:,i][vabs[:,:,i] .> 0.0]), 1:size(vy)[3])
        vabs_error = ((vx_ratio_max .* vx_error).^2 .+ (vy_ratio_max .* vy_error).^2).^0.5
        vabs_error = convert(typeof(vx_error), vabs_error)  
    else         
        vx_error = nothing
        vy_error = nothing
        vabs_error = nothing
    end

    # Run some basic tests 
    nx, ny, ntimes = size(vx)
    @assert length(date1) == length(date2) == ntimes
    @assert nx == ny == 250

    # Spatial preprocessing
    proj_zone = ncgetatt(file, "mapping", "utm_zone_number")
    transform(X,Y) = Sleipnir.UTMercator(X, Y; zone=proj_zone, hemisphere=:north) 
    Coordinates = transform.(x, y)
    latitudes = map(x -> x.lat.val, Coordinates)
    longitudes = map(x -> x.lon.val, Coordinates)

    return SurfaceVelocityData(x=x, y=y, lat=latitudes, lon=longitudes, vx=vx, vy=vy, vabs=vabs, vx_error=vx_error, vy_error=vy_error, vabs_error=vabs_error, date=date_mean, date1=date1, date2=date2, date_error=date_error)
end

"""
    max_or_empyt(A::Array)

Return maximum value for non-empty arrays. 
This is just required to compute the error in the absolute velocity.
"""
function max_or_empty(A::Array)
    if length(A) == 0
        return 0.0
    else
        return maximum(A)
    end
end