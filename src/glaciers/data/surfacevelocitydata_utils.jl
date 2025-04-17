export initialize_surfacevelocitydata

"""
    initialize_surfacevelocitydata(
        raster::Union{String, RasterStack};
        glacier::Union{G, Nothing}=nothing,
        mapping::VM=MeanDateVelocityMapping(),
        compute_vabs_error::Bool=true
    ) where {G <: AbstractGlacier, VM <: VelocityMapping}

Initialize SurfaceVelocityData from Rabatel et. al (2023).

Arguments:
- `raster::Union{String, RasterStack}`: RasterStack or path of the netCDF file with surface velocity data.
- `glacier::Union{G, Nothing}`: Glacier associated to the surface velocity datacube.
    When provided, the surface velocity data are gridded on the glacier grid using
    the `mapping`.
- `mapping::VM`: Mapping to use in order to grid the data from the coordinates of
    the velocity product datacube to the glacier grid.
- `compute_vabs_error::Bool`: Whether to compute the absolute error uncertainty.
"""
function initialize_surfacevelocitydata(
    raster::Union{String, RasterStack};
    glacier::Union{G, Nothing}=nothing,
    mapping::VM=MeanDateVelocityMapping(),
    compute_vabs_error::Bool=true
) where {G <: AbstractGlacier, VM <: VelocityMapping}

    velRast = isa(raster, String) ? RasterStack(raster, lazy=true) : raster

    # Boolean variable indicating if the file represents inporpolated data or not.
    interp = :xcount_x in keys(velRast)

    mapToGlacierGrid = !isnothing(glacier)

    # Date of first adquisition
    date1 = velRast.date1.data[:]
    # Date of second adquisition
    date2 = velRast.date2.data[:]

    # Middle date
    if !interp
        date1 = Int32.(date1)
        date2 = Int32.(date2)
        date_mean = mjd.(0.5 .* date1 .+ 0.5 .* date2)
        date1 = mjd.(date1)
        date2 = mjd.(date2)
    else
        date_mean = dims(velRast, :mid_date).val.data
    end
    date1 = DateTime.(date1)
    date2 = DateTime.(date2)
    date_error = date2 .- date1

    # Read data from netcdf file
    x = dims(velRast, :X).val.data
    y = dims(velRast, :Y).val.data

    # Velocity in the x direction (m/yr)
    vx = velRast.vx.data
    vy = velRast.vy.data

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

    if mapToGlacierGrid
        # Here vx and vy are of type DiskArrays.BroadcastDiskArray
        x, y, vx, vy = grid(glacier, latitudes, longitudes, vx, vy, mapping)
        # Also retrieve the regional coordinates of the glacier grid through x and y
        latitudes = glacier.Coords["lat"]
        longitudes = glacier.Coords["lon"]
        # Here vx and vy have been read from disk and they are arrays
    else
        # Access elements by converting a DiskArrays.BroadcastDiskArray to a real array
        vx = vx[:,:,:]
        vy = vy[:,:,:]
    end

    # Compute absolute velocity
    vabs = (vx.^2 .+ vy.^2).^0.5
    # The sqrt operation in Julia promotes Float32 to Float64. We convert manually
    # to keep type consistency
    vabs = convert(typeof(vx), vabs)

    # The replace function promotes Float32 to Float64. We convert
    # manually to keep consistency
    vx = [eltype(x).(replace(vx[:,:,i], missing => NaN)) for i in 1:size(vx, 3)]
    vy = [eltype(x).(replace(vy[:,:,i], missing => NaN)) for i in 1:size(vy, 3)]
    vabs = [eltype(x).(replace(vabs[:,:,i], missing => NaN)) for i in 1:size(vabs, 3)]

    # Set ice velocity to NaN outside of the glacier outlines
    for i in range(1, size(vx,1))
        vx[i][glacier.H₀ .== 0] .= NaN
        vy[i][glacier.H₀ .== 0] .= NaN
        vabs[i][glacier.H₀ .== 0] .= NaN
    end

    # Error is reported once per timespan, so upper bounds are given by absolute error
    if !interp
        # Since variables below are DiskArrays.BroadcastDiskArray, slice them to access their elements
        vx_error = eltype(x).(replace(velRast.error_vx.data[:], missing => NaN))
        vy_error = eltype(x).(replace(velRast.error_vy.data[:], missing => NaN))
        # Absolute error uncertainty using propagation of uncertanties
        if compute_vabs_error
            vx_ratio_max = map(i -> ratio_max(vx[i], vabs[i]), 1:size(vx,1))
            vy_ratio_max = map(i -> ratio_max(vy[i], vabs[i]), 1:size(vy,1))
            vabs_error = ((vx_ratio_max .* vx_error).^2 .+ (vy_ratio_max .* vy_error).^2).^0.5
            vabs_error = convert(typeof(vx_error[:]), vabs_error[:])
        else
            vabs_error = nothing
        end
    else
        vx_error = nothing
        vy_error = nothing
        vabs_error = nothing
    end

    return SurfaceVelocityData(
        x=x, y=y, lat=latitudes, lon=longitudes,
        vx=vx, vy=vy, vabs=vabs,
        vx_error=vx_error, vy_error=vy_error, vabs_error=vabs_error,
        date=date_mean, date1=date1, date2=date2, date_error=date_error,
        glacierGridded=mapToGlacierGrid
    )
end

"""
    ratio_max(v, vabs)

Compute the maximum ratio between v and vabs at points where the value of vabs is
not a NaN.
"""
function ratio_max(v, vabs)
    mask = replace(v, NaN => 0.0) .> 0.0
    return max_or_empty(abs.(v[mask]) ./ vabs[mask])
end
"""
    max_or_empty(A::Array)

Return maximum value for non-empty arrays.
This is just required to compute the error in the absolute velocity.
"""
function max_or_empty(A::Array)
    return length(A) == 0 ? 0.0 : maximum(A)
end
