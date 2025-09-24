export initialize_surfacevelocitydata, random_spatially_coherent_mask

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
    raster::Union{String, <: RasterStack};
    glacier::Union{G, Nothing}=nothing,
    mapping::VM=MeanDateVelocityMapping(),
    compute_vabs_error::Bool=true
) where {G <: AbstractGlacier, VM <: VelocityMapping}

    velRast = isa(raster, String) ? RasterStack(raster, lazy=true) : raster

    # Boolean variable indicating if the file represents interpolated data or not.
    interp = :xcount_x in keys(velRast)

    mapToGlacierGrid = !isnothing(glacier)

    # Date of first adquisition
    date1 = velRast.date1.data[:]
    # Date of second adquisition
    date2 = velRast.date2.data[:]

    date1 = Int32.(date1)
    date2 = Int32.(date2)

    # Middle date
    if !interp
        date_mean = mjd.(0.5 .* date1 .+ 0.5 .* date2)
        date1 = mjd.(date1)
        date2 = mjd.(date2)
    else
        date_mean = dims(velRast, :mid_date).val.data
        date1 = CFTime.timedecode(date1, metadata(velRast.date1)["units"])
        date2 = CFTime.timedecode(date2, metadata(velRast.date2)["units"])
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

        # Set ice velocity to NaN outside of the glacier outlines
        mask = glacier.H₀ .== 0
        for i in range(1, size(vx,3))
            vx[mask,i] .= missing
            vy[mask,i] .= missing
        end
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
        isGridGlacierAligned=mapToGlacierGrid
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


"""
    grid(
        glacier::G,
        latitudes::Vector{F},
        longitudes::Vector{F},
        vx::Union{FileArray, Array{Union{Missing, F}, 3}},
        vy::Union{FileArray, Array{Union{Missing, F}, 3}},
        mapping::VM
    ) where {
        G <: AbstractGlacier,
        F <: AbstractFloat,
        VM <: VelocityMapping,
        FileArray <: Rasters.FileArray
    }

Grid velocity data onto the glacier grid following the prescribed mapping.
This function maps the 3 dimensional surface velocities (x, y and t) to the glacier grid.
The provided surface velocities can be a `Rasters.FileArray` which happens when the
`RasterStack` is instantiated in lazy mode. In this situation, only the smallest
cube that contains all the needed data to construct the mapping is read from disk.
The returned velocity variables have shape `(nTimes, nx, ny)` where `nTimes` is the
number of time steps and `(nx, ny)` is the size of the glacier grid.

Arguments:
- `glacier::G`: Glacier instance which determines the glacier on which the
    velocities are projected onto.
- `latitudes::Vector{F}`: Vector of latitude values of the original surface
    velocity grid.
- `longitudes::Vector{F}`: Vector of longitude values of the original surface
    velocity grid.
- `vx::Union{FileArray, Array{Union{Missing, F}, 3}}`: X component of the original
    surface velocities. It can be either a `Rasters.FileArray` if the datacube is
    read in lazy mode, or a plain 3 dimensional array.
- `vy::Union{FileArray, Array{Union{Missing, F}, 3}}`: Y component of the original
    surface velocities. It can be either a `Rasters.FileArray` if the datacube is
    read in lazy mode, or a plain 3 dimensional array.
- `mapping::VM`: Mapping to use.

Returns:
- `xG`: A vector that gives the x coordinates of the glacier grid.
- `yG`: A vector that gives the y coordinates of the glacier grid.
- `vxG`: A 3 dimensional array of the x component of the velocity gridded onto the
    glacier grid.
- `vyG`: A 3 dimensional array of the y component of the velocity gridded onto the
    glacier grid.
"""
function grid(
    glacier::G,
    latitudes::Vector{F},
    longitudes::Vector{F},
    vx::Union{FileArray, Array{Union{Missing, F}, 3}},
    vy::Union{FileArray, Array{Union{Missing, F}, 3}},
    mapping::VM
) where {
    G <: AbstractGlacier,
    F <: AbstractFloat,
    VM <: VelocityMapping,
    FileArray <: Rasters.FileArray
}
    params_projection = glacier.params_projection
    transformReverse(lat,lon) = ReverseUTMercator(
        lat, lon;
        k=params_projection["k"],
        cenlon=params_projection["lon_0"], cenlat=params_projection["lat_0"],
        x0=params_projection["x_0"], y0=params_projection["y_0"]
    )
    xG = map(x -> x.x.val, transformReverse.(glacier.cenlat, glacier.Coords["lon"]))
    yG = map(x -> x.y.val, transformReverse.(glacier.Coords["lat"], glacier.cenlon))
    cenlatVel = (latitudes[end]+latitudes[begin])/2
    cenlonVel = (longitudes[end]+longitudes[begin])/2
    xV = map(x -> x.x.val, transformReverse.(cenlatVel, longitudes))
    yV = map(x -> x.y.val, transformReverse.(latitudes, cenlonVel))

    velType = eltype(vx)
    vxG = zeros(velType, glacier.nx, glacier.ny, size(vx,3))
    vyG = zeros(velType, glacier.nx, glacier.ny, size(vx,3))

    if mapping.spatialInterp == :nearest
        # We express each of the coordinates of the velocity grid as (ΔxV*ix+bx, ΔyV*iy+by)
        # While this is an approximation since the coordinates do not truly live on a grid because of the projection, this error does not exceed one meter for most glaciers. This allows us to avoid having to compute an argmin.
        ΔxV = mean(diff(xV))
        ΔyV = mean(diff(yV))
        bx = xV[1]-ΔxV
        by = yV[1]-ΔyV
        ix = collect(1:length(longitudes))
        iy = collect(1:length(latitudes))
        indx = Int.(round.((xG .- bx)./ΔxV))
        indy = Int.(round.((yG .- by)./ΔyV))

        # Lazy arrays need to be read by block, hence we read the smallest block of data that contains all the points we need
        @assert size(vx,1)>=maximum(indx) "It looks like the datacube doesn't cover the whole glacier grid on the x-axis, please check that you use the right datacube for glacier $(glacier.rgi_id)."
        @assert size(vx,2)>=maximum(indy) "It looks like the datacube doesn't cover the whole glacier grid on the y-axis, please check that you use the right datacube for glacier $(glacier.rgi_id)."
        block_vx = vx[minimum(indx):maximum(indx),minimum(indy):maximum(indy),:]
        block_vy = vy[minimum(indx):maximum(indx),minimum(indy):maximum(indy),:]

        # Assign to each point of the glacier grid the closest point on the grid of surface velocities
        shiftx = -minimum(indx)+1
        shifty = -minimum(indy)+1
        for ix in range(1,length(glacier.Coords["lon"])), iy in range(1,length(glacier.Coords["lat"]))
            vxG[ix,iy,begin:end] .= block_vx[indx[ix]+shiftx,indy[iy]+shifty,:]
            vyG[ix,iy,begin:end] .= block_vy[indx[ix]+shiftx,indy[iy]+shifty,:]
        end
    else
        throw("$(mapping.spatialInterp) spatial interpolation method is not implemented")
    end
    return xG, yG, vxG, vyG
end

"""
    random_spatially_coherent_mask(h::Integer, w::Integer; sigma::Real=1.0, threshold::Real=0.0) -> BitMatrix
    random_spatially_coherent_mask(mask::BitMatrix; sigma::Real=1.0, threshold::Real=0.0) -> BitMatrix

Generate a random binary mask with **spatially correlated patches** rather than pixel-wise independent noise.
This is done by drawing white noise, applying a Gaussian low-pass filter in the frequency domain, and thresholding
the result.

# Arguments
- `h::Integer`, `w::Integer`: Height and width of the mask.
- `mask::BitMatrix`: An existing binary mask. The generated spatially coherent mask
    will be applied elementwise (`.&`) to this mask.
- `sigma::Real=1.0`: Controls the spatial correlation length. Larger values produce
    smoother, larger patches.
- `threshold::Real=0.0`: Threshold applied to the filtered noise. Higher values
    result in sparser masks. Statistically, setting the threshold to zero results in
    a mask with half pixels to true.

# Returns
A `BitMatrix` of size `(h, w)` containing `true` in patchy regions and `false` elsewhere.

# Examples
```julia
# Generate a new 256×256 patchy mask
mask = random_spatially_coherent_mask(256, 256; sigma=8.0, threshold=0.0)

# Apply patchy masking to an existing mask
base = trues(128, 128)
patchy = random_spatially_coherent_mask(base; sigma=5.0, threshold=0.3)
"""
function random_spatially_coherent_mask(h::Integer, w::Integer; sigma::Real=1.0, threshold::Real=0.0)
    # 1) white noise
    noise = randn(h, w)

    # 2) forward 2D rfft
    f = fft(noise)   # size (h, w÷2 + 1)

    # 3) build 2D frequency grid
    kx = fftfreq(h)                  # length h
    ky = fftfreq(w) # (0:(w÷2)) ./ w              # length w÷2 + 1 (since rfft2 keeps nonnegative y-freqs)

    KX = reshape(kx, h, 1) .* ones(1, length(ky))  # h × (w÷2+1)
    KY = ones(h, 1) .* reshape(ky, 1, length(ky))  # h × (w÷2+1)

    # 4) Gaussian low-pass filter
    power = exp.(- (KX.^2 .+ KY.^2) .* (2π*sigma)^2)

    # 5) apply filter in frequency space
    f_filtered = f .* power

    # 6) invert back to real space
    smooth = real(ifft(f_filtered))  # explicitly give original size

    # 7) threshold for binary mask
    return smooth .> threshold
end
function random_spatially_coherent_mask(mask::BitMatrix; sigma::Real=1.0, threshold::Real=0.0)
    h, w = size(mask)
    random_mask = random_spatially_coherent_mask(h, w; sigma=sigma, threshold=threshold)
    return mask .& random_mask
end
