export initialize_surfacevelocitydata

"""
    initialize_surfacevelocitydata(file::String; interp::Bool=false, compute_vabs_error::Bool=true)

Initialize SurfaceVelocityData from Rabatel et. al (2023).

Arguments:
- `file::String`: Path of the netCDF file with surface velocity data.
- `interp::Bool`: Boolean variable indicating if the file represents inporpolated data or not.
- `compute_vabs_error::Bool`: Whether to compute the absolute error uncertainty.
"""
function initialize_surfacevelocitydata(file::String; interp::Bool=false, compute_vabs_error::Bool=true)

    velRast = RasterStack(file)

    # Date of first adquisition
    date1 = velRast.date1.data
    # Date of second adquisition
    date2 = velRast.date2.data

    # Middle date
    if !interp
        # Convert to Julian datetime
        date1_julian = (datetime2julian.(date1) .- 2400000.5)
        date2_julian = (datetime2julian.(date2) .- 2400000.5)
        # Convert dates from Modified Julian Days to datetimes
        date_mean = mjd.(0.5 .* date1_julian .+ 0.5 .* date2_julian)
    else
        date_mean = dims(velRast, :mid_date).val.data
    end
    date_error = date2 .- date1

    # Read data from netcdf file
    x = dims(velRast, :X).val.data
    y = dims(velRast, :Y).val.data

    # Velocity in the x direction (m/yr)
    vx = velRast.vx.data
    vy = velRast.vy.data

    # Compute absolute velocity
    vabs = (vx.^2 .+ vy.^2).^0.5
    # The sqrt operation in Julia promotes Float32 to Float64. We convert manually
    # to keep type consistency
    vabs = convert(typeof(vx), vabs)

    # Error is reported once per timespan, so upper bounds are given by absolute error
    if !interp
        vx_error = velRast.error_vx.data
        vy_error = velRast.error_vy.data
        # Absolute error uncertainty using propagation of uncertanties
        if compute_vabs_error
            vx_ratio_max = map(i -> max_or_empty(abs.(vx[:,:,i][vabs[:,:,i] .> 0.0]) ./ vabs[:,:,i][vabs[:,:,i] .> 0.0]), 1:size(vx)[3])
            vy_ratio_max = map(i -> max_or_empty(abs.(vy[:,:,i][vabs[:,:,i] .> 0.0]) ./ vabs[:,:,i][vabs[:,:,i] .> 0.0]), 1:size(vy)[3])
            vabs_error = ((vx_ratio_max .* vx_error).^2 .+ (vy_ratio_max .* vy_error).^2).^0.5
            vabs_error = convert(typeof(vx_error), vabs_error)
        else
            vabs_error = nothing
        end
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
    params_projection = parse_proj(metadata(velRast)["proj4"])
    hemisphere = (y[end]+y[begin])/2 >= 0 ? :north : :south
    transform(X,Y) = Sleipnir.UTMercator(X, Y; zone=Int(params_projection["zone"]), hemisphere=hemisphere)
    latitudes = map(x -> x.lat.val, transform.(mean(x), y))
    longitudes = map(x -> x.lon.val, transform.(x, mean(y)))

    return SurfaceVelocityData(x=x, y=y, lat=latitudes, lon=longitudes, vx=vx, vy=vy, vabs=vabs, vx_error=vx_error, vy_error=vy_error, vabs_error=vabs_error, date=date_mean, date1=date1, date2=date2, date_error=date_error)
end

"""
    max_or_empty(A::Array)

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
