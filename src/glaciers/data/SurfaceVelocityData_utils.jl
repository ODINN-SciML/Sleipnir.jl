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
    compute_vabs_error::Bool=true,
    flag::Union{String, <: RasterStack, Nothing} = nothing
) where {G <: AbstractGlacier, VM <: VelocityMapping}

    velRast = isa(raster, String) ? RasterStack(raster, lazy=true) : raster

    # Boolean variable indicating if the file represents interpolated data or not.
    interp = :xcount_x in keys(velRast)

    mapToGlacierGrid = !isnothing(glacier)

    # Dates in the dataset are integers representing a single day (no hour or minutes information)
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

    # Velocity in the x / longiudinal direction (m/yr)
    # Positive values correspond to W-E orientation
    vx = velRast.vx.data
    # Velocity in the y / latitudinal direction (m/yr)
    # Positive values correspond to S-N orientation
    vy = velRast.vy.data

    # Run some basic tests
    nx, ny, ntimes = size(vx)
    @assert length(date1) == length(date2) == ntimes
    @assert nx == ny == 250

    # Spatial preprocessing
    params_projection = parse_proj(metadata(velRast)["proj4"])
    hemisphere = (y[end]+y[begin])/2 >= 0 ? :north : :south
    transform(X,Y) = Sleipnir.UTMercator(X, Y; zone=Int(params_projection["zone"]), hemisphere=hemisphere)
    latitudes_vel = map(x -> x.lat.val, transform.(mean(x), y))
    longitudes_vel = map(x -> x.lon.val, transform.(x, mean(y)))

    # Ice velocity mask product
    vflag = isnothing(flag) ? nothing : initialize_surfacevelocitydata_mask(velRast, flag)

    if mapToGlacierGrid
        # Also retrieve the regional coordinates of the glacier grid through x and y
        latitudes_glacier = glacier.Coords["lat"]
        longitudes_glacier = glacier.Coords["lon"]

        # Define mask where there is velocity data
        min_distance = 0.5 * (glacier.Δx^2 + glacier.Δy^2)^0.5
        Latitudes_vel = latitudes_vel * ones(length(longitudes_vel))'
        Longitudes_vel = ones(length(latitudes_vel)) * longitudes_vel'
        mask_data = [
            minimum(local_distance.(Latitudes_vel, Longitudes_vel, Ref(lat), Ref(lon))) .> min_distance
            for lon in longitudes_glacier,
            lat in latitudes_glacier
            ]
        # Define mask where there is no ice
        mask_ice = glacier.H₀ .== 0

        # Test that at least part of the target glacier falls inside datacube
        # velocity datacube corners
        lat_v_begin, lat_v_end = minimum(latitudes_vel), maximum(latitudes_vel)
        lon_v_begin, lon_v_end = minimum(longitudes_vel), maximum(longitudes_vel)
        # glacier corners
        is_icy = glacier.H₀ .> 0.0
        icy_latitudes_glacier = latitudes_glacier[any.(eachcol(is_icy))]
        icy_longitudes_glacier = longitudes_glacier[any.(eachrow(is_icy))]

        lat_g_begin, lat_g_end = minimum(icy_latitudes_glacier), maximum(icy_latitudes_glacier)
        lon_g_begin, lon_g_end = minimum(icy_longitudes_glacier), maximum(icy_longitudes_glacier)

        # Overlapping occurs when one of the glacier corners is inside velocity datacube:
        lw, up = (lat_v_begin, lon_v_begin), (lat_v_end, lon_v_end)
        glacier_in_datacube = [
            all(lw .<= [lat_g_begin, lon_g_begin] .<= up),
            all(lw .<= [lat_g_begin, lon_g_end] .<= up),
            all(lw .<= [lat_g_end, lon_g_begin] .<= up),
            all(lw .<= [lat_g_end, lon_g_end] .<= up)
        ]
        @assert any(glacier_in_datacube) "Datacube doesn't include any region of the glacier, please check that you use the right datacube for glacier $(glacier.rgi_id)."
        if any(.!glacier_in_datacube)
            coverage = round(100 * sum(.!mask_ice .&& .!mask_data) / sum(.!mask_ice); digits = 2)
            @assert coverage > 0.0 "Datacube doesn't include any region of the glacier, please check that you use the right datacube for glacier $(glacier.rgi_id)."
            @warn "Glacier is not enterely included in the datacube. Current datacube covers $(coverage)% of glacier $(glacier.rgi_id)."
        end

        # Here vx and vy are of type DiskArrays.BroadcastDiskArray
        x, y, vx, vy, vflag = grid(glacier, latitudes_vel, longitudes_vel, vx, vy, vflag, mapping)
        # Set ice velocity to NaN outside of the glacier outlines
        mask = mask_ice .|| mask_data
        for i in range(1, size(vx,3))
            vx[mask, i] .= missing
            vy[mask, i] .= missing
        end
    else
        # Access elements by converting a DiskArrays.BroadcastDiskArray to a real array
        vx = vx[:, :, :]
        vy = vy[:, :, :]
    end

    # Define coordinates used for the velocity dataset
    latitudes = mapToGlacierGrid ? latitudes_glacier : latitudes_vel
    longitudes = mapToGlacierGrid ? longitudes_glacier : longitudes_vel

    # Align ice surface velocity vector with glacier convention of easting and northing
    # Glacier convention follows matrix index notation: first index is latitude (northing), second index is longitude (easting)
    # The ice surface velocity product reports velocities vx and vy in east-west and south-north orientations, respectivelly
    # To check for the right orientation of the velocity vectors, we see how coordinates are encoded in each data object
    vx = all(diff(latitudes) .> 0.0) ? vx : .- vx
    vy = all(diff(longitudes) .> 0.0) ? vy : .- vy

    # Compute absolute velocity
    vabs = (vx.^2 .+ vy.^2).^0.5
    # The sqrt operation in Julia promotes Float32 to Float64. We convert manually
    # to keep type consistency
    vabs = convert(typeof(vx), vabs)

    # The replace function promotes Float32 to Float64. We convert
    # manually to keep consistency
    vx = [eltype(x).(replace(vx[:, :, i], missing => NaN)) for i in 1:size(vx, 3)]
    vy = [eltype(x).(replace(vy[:, :, i], missing => NaN)) for i in 1:size(vy, 3)]
    vabs = [eltype(x).(replace(vabs[:, :, i], missing => NaN)) for i in 1:size(vabs, 3)]

    # Error is reported once per timespan, so upper bounds are given by absolute error
    if !interp
        # Since variables below are DiskArrays.BroadcastDiskArray, slice them to access their elements
        vx_error = eltype(x).(replace(velRast.error_vx.data[:], missing => NaN))
        vy_error = eltype(x).(replace(velRast.error_vy.data[:], missing => NaN))
        # Absolute error uncertainty using propagation of uncertanties
        if compute_vabs_error
            vx_ratio_max = map(i -> ratio_max(vx[i], vabs[i]), 1:size(vx, 1))
            vy_ratio_max = map(i -> ratio_max(vy[i], vabs[i]), 1:size(vy, 1))
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
        x = x, y = y,
        lat = latitudes, lon = longitudes,
        vx = vx, vy = vy,
        vabs = vabs,
        vx_error = vx_error, vy_error = vy_error, vabs_error = vabs_error,
        date = date_mean, date1 = date1, date2 = date2, date_error = date_error,
        flag = vflag,
        isGridGlacierAligned = mapToGlacierGrid
    )
end

function initialize_surfacevelocitydata_mask(
    data::RasterStack,
    flag::Union{String, <: RasterStack, Nothing} = nothing
)
    flagRast = isa(flag, String) ? RasterStack(flag, lazy = true) : flag

    # Subset mask based on glacier datacube
    # Flag is centered in the middle pixel, so we shift It
    Δ = 26.0
    Xs = dims(data, :X).val.data
    X_begin, X_end = minimum(Xs) - Δ, maximum(Xs) + Δ
    Ys = dims(data, :Y).val.data
    Y_begin, Y_end = minimum(Ys) - Δ, maximum(Ys) + Δ

    flagRast_subset = flagRast[X(X_begin..X_end), Y(Y_begin..Y_end)]

    # We align the axes of the glacier subset with the flag
    X_glacier_increasing = first(Xs) < last(Xs)
    Y_glacier_increasing = first(Ys) < last(Ys)
    X_flag_increasing = first(dims(flagRast_subset, :X).val.data) < last(dims(flagRast_subset, :X).val.data)
    Y_flag_increasing = first(dims(flagRast_subset, :Y).val.data) < last(dims(flagRast_subset, :Y).val.data) 

    X_reverse = xor(X_glacier_increasing, X_flag_increasing)
    Y_reverse = xor(Y_glacier_increasing, Y_flag_increasing)
    flagRast_subset = X_reverse ? reverse(flagRast_subset, dims = X) : flagRast_subset
    flagRast_subset = Y_reverse ? reverse(flagRast_subset, dims = Y) : flagRast_subset

    # Unreliable pixels are indicated with 0
    # Reliable pixels are indicared with 1
    return flagRast_subset.layer1.data .== 0
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
    vflag::Union{BitMatrix, Nothing},
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
    # Glacier coordinates in northing/easting coordinates
    xG = map(x -> x.x.val, transformReverse.(glacier.cenlat, glacier.Coords["lon"]))
    yG = map(x -> x.y.val, transformReverse.(glacier.Coords["lat"], glacier.cenlon))

    # Velocity coordinates in northing/easting coordinates
    cenlatVel = (latitudes[end]+latitudes[begin])/2
    cenlonVel = (longitudes[end]+longitudes[begin])/2
    xV = map(x -> x.x.val, transformReverse.(cenlatVel, longitudes))
    yV = map(x -> x.y.val, transformReverse.(latitudes, cenlonVel))

    # Velocity tensor in glacier grid
    velType = eltype(vx)
    vxG = zeros(velType, glacier.nx, glacier.ny, size(vx,3))
    vyG = zeros(velType, glacier.nx, glacier.ny, size(vx,3))
    flagG = isnothing(vflag) ? nothing : BitMatrix(falses(glacier.nx, glacier.ny))

    if mapping.spatialInterp == :nearest
        # We express each of the coordinates of the velocity grid as (ΔxV*ix+bx, ΔyV*iy+by)
        # While this is an approximation since the coordinates do not truly live on a grid because of the projection, this error does not exceed one meter for most glaciers. This allows us to avoid having to compute an argmin.
        ΔxV = mean(diff(xV))
        ΔyV = mean(diff(yV))
        bx = xV[1] - ΔxV
        by = yV[1] - ΔyV
        ix = collect(1:length(longitudes))
        iy = collect(1:length(latitudes))
        indx = Int.(round.((xG .- bx)./ΔxV))
        indy = Int.(round.((yG .- by)./ΔyV))

        # Lazy arrays need to be read by block, hence we read the smallest block of data that contains all the points we need
        indx_lw = max(minimum(indx), 1)
        indx_up = min(maximum(indx), size(vx, 1))
        indy_lw = max(minimum(indy), 1)
        indy_up = min(maximum(indy), size(vx, 1))
        block_vx = vx[indx_lw:indx_up, indy_lw:indy_up, :]
        block_vy = vy[indx_lw:indx_up, indy_lw:indy_up, :]
        block_flag = isnothing(vflag) ? nothing : vflag[indx_lw:indx_up, indy_lw:indy_up]

        # Assign to each point of the glacier grid the closest point on the grid of surface velocities
        shiftx = - indx_lw + 1
        shifty = - indy_lw + 1
        for ix in range(1, length(glacier.Coords["lon"])), iy in range(1, length(glacier.Coords["lat"]))
            ixv, iyv = indx[ix] + shiftx, indy[iy] + shifty
            if checkbounds(Bool, block_vx, ixv, iyv, 1)
                vxG[ix,iy,begin:end] .= block_vx[ixv, iyv, :]
                vyG[ix,iy,begin:end] .= block_vy[ixv, iyv, :]
                if !isnothing(vflag)
                    flagG[ix, iy] = block_flag[ixv, iyv]
                end
            end
        end
    else
        throw("$(mapping.spatialInterp) spatial interpolation method is not implemented")
    end
    return xG, yG, vxG, vyG, flagG
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
"""
Compute distance between coordinates in meters
"""
function local_distance(lat1, lon1, lat2, lon2)
    # average latitude for scaling longitude
    φ = deg2rad((lat1 + lat2) / 2)
    dx = (lon2 - lon1) * 111320 * cos(φ)  # meters in x (east-west)
    dy = (lat2 - lat1) * 111320           # meters in y (north-south)
    return sqrt(dx^2 + dy^2)
end
"""
    combine_velocity_data(refVelocities; merge=false)

Combine multiple ice surface velocity datasets into a single `SurfaceVelocityData` object.

# Arguments
- `refVelocities::Vector{SurfaceVelocityData}`: A vector of ice surface velocity datasets to combine. Each element must have the same grid alignment (`isGridGlacierAligned` must be `true` for all).
- `merge::Bool=false`: If `true`, velocities with the same `date` are averaged, and corresponding date ranges (`date1`, `date2`) are reduced to their min/max. If `false`, data are simply concatenated.

# Returns
- `SurfaceVelocityData`: A single object containing the combined velocity data, including `vx`, `vy`, `vabs` and their associated errors, as well as coordinate (`x`, `y`, `lat`, `lon`) and date information. The `isGridGlacierAligned` field reflects whether all input datasets were aligned.

# Notes
- The function asserts that all input datasets are aligned on the same grid.
- When `merge=true`, velocities and errors are averaged over datasets sharing the same `date`.
- Uses `nanmean` for averaging to handle missing data.
- `date_error` is set to `nothing` when merging.
"""
function combine_velocity_data(refVelocities; merge = false)
    # Check all surfaces are on the same grid as the glacier
    isGridGlacierAligned = [refV.isGridGlacierAligned for refV in refVelocities]
    @assert all(isGridGlacierAligned) "Different ice surface velocity datasets are not alligned"

    # Shared features
    x = refVelocities[begin].x
    y = refVelocities[begin].y
    lat = refVelocities[begin].lat
    lon = refVelocities[begin].lon
    # Concat features
    vx = reduce(vcat, [refV.vx for refV in refVelocities])
    vy = reduce(vcat, [refV.vy for refV in refVelocities])
    vabs = reduce(vcat, [refV.vabs for refV in refVelocities])
    vx_error = reduce(vcat, [refV.vx_error for refV in refVelocities])
    vy_error = reduce(vcat, [refV.vy_error for refV in refVelocities])
    vabs_error = reduce(vcat, [refV.vabs_error for refV in refVelocities])
    date = reduce(vcat, [refV.date for refV in refVelocities])
    date1 = reduce(vcat, [refV.date1 for refV in refVelocities])
    date2 = reduce(vcat, [refV.date2 for refV in refVelocities])
    date_error = reduce(vcat, [refV.date_error for refV in refVelocities])
    # We defined the combined mask based if the flag is activated in one of the datacubes
    flag = reduce(.|, [refV.flag for refV in refVelocities])

    if merge
        date_unique = sort(unique(date))
        date1_unique = similar(date1, length(date_unique))
        date2_unique = similar(date2, length(date_unique))
        vx_unique = similar(vx, length(date_unique))
        vy_unique = similar(vy, length(date_unique))
        vabs_unique = similar(vabs, length(date_unique))
        vx_error_unique = similar(vx_error, length(date_unique))
        vy_error_unique = similar(vy_error, length(date_unique))
        vabs_error_unique = similar(vabs_error, length(date_unique))

        for (i, dt) in enumerate(date_unique)
            # Reduce datasets for unique datetime
            vx_unique[i] = mapslices(
                nanmean,
                cat(vx[date .== dt]..., dims = 3);
                dims = 3
                )[:, :, 1]
            vy_unique[i] = mapslices(
                nanmean,
                cat(vy[date .== dt]..., dims = 3);
                dims = 3
                )[:, :, 1]
            vabs_unique[i] = mapslices(
                nanmean,
                cat(vabs[date .== dt]..., dims = 3);
                dims = 3
                )[:, :, 1]
            vx_error_unique[i] = mapslices(
                nanmean,
                cat(vx_error[date .== dt]..., dims = 3);
                dims = 3
                ) |> only
            vy_error_unique[i] = mapslices(
                nanmean,
                cat(vy_error[date .== dt]..., dims = 3);
                dims = 3
                ) |> only
            vabs_error_unique[i] = mapslices(
                nanmean,
                cat(vabs_error[date .== dt]..., dims = 3);
                dims = 3
                ) |> only
            date1_unique[i] = mapslices(
                minimum,
                cat(date1[date .== dt]..., dims = 3);
                dims = 3
                ) |> only
            date2_unique[i] = mapslices(
                maximum,
                cat(date2[date .== dt]..., dims = 3);
                dims = 3
                ) |> only
        end
        return SurfaceVelocityData(
            x = x,
            y = y,
            lat = lat,
            lon = lon,
            vx = vx_unique,
            vy = vy_unique,
            vabs = vabs_unique,
            vx_error = vx_error_unique,
            vy_error = vy_error_unique,
            vabs_error = vabs_error_unique,
            date = date_unique,
            date1 = date1_unique,
            date2 = date2_unique,
            date_error = nothing,
            flag = flag,
            isGridGlacierAligned = all(isGridGlacierAligned)
        )
    else
        return SurfaceVelocityData(
            x = x,
            y = y,
            lat = lat,
            lon = lon,
            vx = vx,
            vy = vy,
            vabs = vabs,
            vx_error = vx_error,
            vy_error = vy_error,
            vabs_error = vabs_error,
            date = date,
            date1 = date1,
            date2 = date2,
            date_error = date_error,
            flag = flag,
            isGridGlacierAligned = all(isGridGlacierAligned)
        )
    end
end
