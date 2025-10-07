
export initialize_glaciers
export is_in_glacier
export glacierName

###############################################
############  FUNCTIONS   #####################
###############################################
"""
    initialize_glaciers(
        rgi_ids::Vector{String},
        params::Parameters;
        velocityDatacubes::Union{Dict{String, String}, Dict{String, RasterStack}}=Dict(),
    )

Initialize glaciers based on provided RGI IDs and parameters.

# Arguments
- `rgi_ids::Vector{String}`: A vector of RGI IDs representing the glaciers to be initialized.
- `params::Parameters`: A `Parameters` object containing simulation parameters.
- `velocityDatacubes::Union{Dict{String, String}, Dict{String, RasterStack}}`: A dictionary that provides for each RGI ID either the path to the datacube or the `RasterStack` with velocity data.

# Returns
- `glaciers::Vector{Glacier2D}`: A vector of initialized `Glacier2D` objects.

# Description
This function performs the following steps:
1. Generates a file for missing glaciers if it does not already exist.
2. Filters out missing glaciers from the provided RGI IDs.
3. Generates raw climate data for the glaciers if necessary.
4. Initializes the glaciers using the provided RGI IDs and parameters.
5. If `use_glathida_data` is enabled in the simulation parameters, assigns GlaThiDa data to the glaciers.

# Errors
- Throws an error if none of the provided RGI IDs have GlaThiDa data.

# Warnings
- Issues a warning if not all glaciers have GlaThiDa data available.

# Example
```julia
# We declare a list of glaciers to be initialized with their RGI IDs
rgi_ids = ["RGI60-11.03638", "RGI60-11.01450", "RGI60-11.02346", "RGI60-08.00203"]
# We initialize those glaciers based on the RGI IDs and the parameters we previously specified
glaciers = initialize_glaciers(rgi_ids, params)
```
"""
function initialize_glaciers(
    rgi_ids::Vector{String},
    params::Parameters;
    velocityDatacubes::Union{Dict{String, String}, Dict{String, Vector{String}}, Dict{String, <: RasterStack}}=Dict{String,String}(),
)

    # Generate missing glaciers file
    missing_glaciers_path = joinpath(params.simulation.working_dir, "data")
    if !isdir(missing_glaciers_path)
        mkdir(missing_glaciers_path)
    end

    if !isfile(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))
        missing_glaciers = Vector([])
        jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"); missing_glaciers)
    end
    filter_missing_glaciers!(rgi_ids, params)

    # Generate raw climate data if necessary
    if params.simulation.test_mode
        # Avoid GitHub CI issue
        map((rgi_id) -> generate_raw_climate_files(rgi_id, params.simulation), rgi_ids)
    else
        pmap((rgi_id) -> generate_raw_climate_files(rgi_id, params.simulation), rgi_ids)
    end

    glaciers = pmap(
        (rgi_id) -> initialize_glacier(rgi_id, params; velocityDatacubes=velocityDatacubes),
        rgi_ids
    )

    if params.simulation.use_glathida_data

        # Obtain H_glathida values for the valid RGI IDs
        H_glathida_values, valid_glaciers = get_glathida!(glaciers, params)
        valid_rgi_ids = [glacier.rgi_id for glacier in valid_glaciers]

        if isempty(valid_rgi_ids)
            error("None of the provided RGI IDs have GlaThiDa.")
        end

        if length(valid_rgi_ids) < length(rgi_ids)
            @warn "Not all glaciers have GlaThiDa data available."
        end

        # Create a mapping from RGI ID to H_glathida value
        rgi_to_H_glathida = Dict(zip(valid_rgi_ids, H_glathida_values))

        # Assign H_glathida to glaciers with valid RGI IDs
        for glacier in glaciers
            if glacier.rgi_id in valid_rgi_ids
                glacier.H_glathida = rgi_to_H_glathida[glacier.rgi_id]
            end
        end
    end

    return glaciers
end

"""
    initialize_glacier(rgi_id::String, parameters::Parameters; smoothing=false)

Initialize a glacier with the given RGI ID and parameters.

# Arguments
- `rgi_id::String`: The RGI (Randolph Glacier Inventory) ID of the glacier.
- `parameters::Parameters`: A struct containing various parameters required for initialization.
- `smoothing::Bool`: Optional. If `true`, apply smoothing to the initial topography. Default is `false`.
- `masking::Union{Int, Nothing, Matrix}`: Type of mask applied to the glacier to determine regions with no ice.
- `velocityDatacubes::Union{Dict{String, String}, Dict{String, RasterStack}}`: A dictionary that provides for each RGI ID either the path to the datacube or the `RasterStack` with velocity data.

# Returns
- `glacier`: An initialized glacier object containing the initial topography and climate data.
"""
function initialize_glacier(
    rgi_id::String,
    parameters::Parameters;
    smoothing::Bool = false,
    masking::Union{Int, Nothing, Matrix} = 2,
    velocityDatacubes::Union{Dict{String, String}, Dict{String, Vector{String}}, Dict{String, <: RasterStack}} = Dict{String,String}(),
)
    # Build glacier and its associated climate
    glacier = Glacier2D(rgi_id, parameters; masking = masking, smoothing = smoothing)

    if get(velocityDatacubes, glacier.rgi_id, "") != ""
        mapping = parameters.simulation.mapping
        datacubes = velocityDatacubes[glacier.rgi_id]
        refVelocities = []
        for datacube in datacubes
            _refVelocity = initialize_surfacevelocitydata(
                datacube;
                glacier = glacier,
                mapping = mapping
                )
            push!(refVelocities, _refVelocity)
        end
        # Rebuild glacier since we cannot change type of `glacier.velocityData`
        # Combine velocity data with unique datetime
        refVelocity = combine_velocity_data(refVelocities; merge = true)
        glacier = Glacier2D(
            glacier,
            velocityData = refVelocity
            )
    end

    return glacier
end

function convertRasterStackToFloat64(rs::RasterStack)
    layerNames = names(rs)
    return RasterStack(
        NamedTuple{Tuple(layerNames)}([Float64.(rs[n]) for n in layerNames]),
        metadata=metadata(rs)
    )
end

"""
    Glacier2D(rgi_id::String, params::Parameters; smoothing=false, test=false)

Build glacier object for a given RGI ID and parameters.

# Arguments
- `rgi_id::String`: The RGI ID of the glacier.
- `params::Parameters`: A `Parameters` object containing simulation parameters.
- `masking::Union{Int, Nothing, Matrix}`: Type of mask applied to the glacier to determine regions with no ice.
- `smoothing::Bool=false`: Optional; whether to apply smoothing to the initial ice thickness. Default is `false`.
- `test::Bool=false`: Optional; test flag. Default is `false`.

# Returns
- `glacier::Glacier2D`: A `Glacier2D` object initialized with the glacier data.

# Description
This function loads and initializes the glacier data for a given RGI ID. It retrieves the initial ice thickness conditions based on the specified source in the parameters, applies optional smoothing, and initializes the glacier's topographical and velocity data. The function also handles Mercator projection for the glacier coordinates and filters glacier borders in high elevations to avoid overflow problems.

# Notes
- The function reverses the matrices for ice thickness, bedrock, and other data to match the required orientation.
- If the Mercator projection includes latitudes larger than 80°, a warning is issued.
- If the glacier data is missing, the function updates a list of missing glaciers and issues a warning.
"""
function Glacier2D(
    rgi_id::String,
    params::Parameters;
    masking::Union{Int, Nothing, BitMatrix} = 2,
    smoothing=false
    )

    # Load glacier gridded data
    F = Sleipnir.Float
    rgi_path = joinpath(prepro_dir, params.simulation.rgi_paths[rgi_id])
    glacier_gd = RasterStack(joinpath(rgi_path, "gridded_data.nc"))
    if Sleipnir.doublePrec
        glacier_gd = convertRasterStackToFloat64(glacier_gd)
    end
    glacier_grid = JSON.parsefile(joinpath(rgi_path, "glacier_grid.json"))
    # Retrieve initial conditions from OGGM
    # initial ice thickness conditions for forward model
    if params.simulation.ice_thickness_source == "Millan22" && params.simulation.use_velocities
        H₀ = F.(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_ice_thickness.data, 0.0))
    elseif params.simulation.ice_thickness_source == "Farinotti19"
        H₀ = F.(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.consensus_ice_thickness.data, 0.0))
    end
    fillNaN!(H₀) # Fill NaNs with 0s to have real boundary conditions
    if smoothing
        println("Smoothing is being applied to initial condition.")
        smooth!(H₀)  # Smooth initial ice thickness to help the solver
    end
    if params.simulation.gridScalingFactor > 1
        H₀ = block_average_pad_edge(H₀, params.simulation.gridScalingFactor)
    end

    try
        # We filter glacier borders in high elevations to avoid overflow problems
        dist_border::Matrix{Sleipnir.Float} = glacier_gd.dis_from_border.data
        if params.simulation.gridScalingFactor > 1
            # Note: this is not mathematically correct and this should be fixed in the future, however since this option is used only in the tests it isn't critical
            dist_border = block_average_pad_edge(dist_border, params.simulation.gridScalingFactor)
        end

        # Define mask where ice can exist (H > 0)
        mask = @match masking begin
            ::Nothing => trues(size(H₀)...)
            ::Int => is_in_glacier(H₀, -masking)
            ::BitMatrix => masking
        end

        nx, ny = params.simulation.gridScalingFactor > 1 ? size(H₀) : glacier_grid["nxny"]

        # Mercator Projection
        params_projection::Dict{String, Float64} = parse_proj(glacier_grid["proj"])
        transform(X,Y) = UTMercator(
            X, Y;
            k=params_projection["k"],
            cenlon=params_projection["lon_0"], cenlat=params_projection["lat_0"],
            x0=params_projection["x_0"], y0=params_projection["y_0"]
        )
        easting = dims(glacier_gd, 1).val
        northing = dims(glacier_gd, 2).val
        latitudes = map(x -> x.lat.val, transform.(Ref(mean(easting)), northing))
        longitudes = map(x -> x.lon.val, transform.(easting, Ref(mean(northing))))
        cenlon::Sleipnir.Float = longitudes[Int(round(nx/2))]
        cenlat::Sleipnir.Float = latitudes[Int(round(ny/2))]
        if maximum(abs.(latitudes)) > 80
            @warn "Mercator projection can fail in high-latitude regions. You glacier includes latitudes larger than 80°."
        end

        S::Matrix{Sleipnir.Float} = glacier_gd.topo.data
        if params.simulation.gridScalingFactor > 1
            S = block_average_pad_edge(S, params.simulation.gridScalingFactor)
            longitudes = longitudes[begin:params.simulation.gridScalingFactor:end]
            latitudes = latitudes[begin:params.simulation.gridScalingFactor:end]
        end
        B = S .- H₀ # bedrock

        Coords = Dict{String,Vector{Float64}}("lon"=> longitudes, "lat"=> latitudes)

        if params.simulation.use_velocities
            V = ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_v.data, 0.0)
            Vx = ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_vx.data, 0.0)
            Vy = ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_vy.data, 0.0)
            fillNaN!(V)
            fillNaN!(Vx)
            fillNaN!(Vy)
        else
            V = zeros(F, size(H₀))
            Vx = zeros(F, size(H₀))
            Vy = zeros(F, size(H₀))
        end
        Δx::Sleipnir.Float = abs.(glacier_grid["dxdy"][1])
        Δy::Sleipnir.Float = abs.(glacier_grid["dxdy"][2])
        if params.simulation.gridScalingFactor > 1
            Δx *= params.simulation.gridScalingFactor
            Δy *= params.simulation.gridScalingFactor
        end
        # Local slope based on smoothed topography
        slope::Matrix{Sleipnir.Float} = glacier_gd.slope.data
        name = get(get_rgi_names(), rgi_id, "")

        # Initialize glacier climate
        climate = Climate2D(rgi_id, params, S, Coords)

        return Glacier2D(
            rgi_id = rgi_id,
            name = name,
            climate = climate,
            H₀ = H₀, S = S, B = B, V = V, Vx = Vx, Vy = Vy,
            A = Sleipnir.Float(4e-17), C = Sleipnir.Float(0.0), n = Sleipnir.Float(3.0),
            slope = slope, dist_border = dist_border, mask = mask,
            Coords = Coords, Δx = Δx, Δy = Δy, nx = nx, ny = ny,
            cenlon = cenlon, cenlat = cenlat,
            params_projection = params_projection
        )

    catch error
        @show error
        missing_glaciers = load(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
        push!(missing_glaciers, rgi_id)
        jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"); missing_glaciers)
        @warn "Glacier without data: $rgi_id. Updating list of missing glaciers. Please try again."
    end
end

# [Begin] Glathida Utilities
"""
    get_glathida!(glaciers::Vector{G}, params::Parameters; force=false) where {G <: Glacier2D}

Retrieve and process glacier thickness data for a vector of `Glacier2D` objects.

# Arguments
- `glaciers::Vector{Glacier2D}`: A vector of `Glacier2D` objects for which the glacier thickness data is to be retrieved.
- `params::Parameters`: A `Parameters` object containing simulation parameters.
- `force::Bool=false`: A boolean flag indicating whether to force the retrieval of glacier thickness data.

# Returns
- `gtd_grids::Vector`: A vector of glacier thickness data grids.
- `glaciers::Vector{Glacier2D}`: The updated vector of `Glacier2D` objects after removing glaciers with no data.

# Description
This function retrieves glacier thickness data for each glacier in the input vector using parallel processing. It updates a list of missing glaciers if any glacier has all data points equal to zero. The function then removes glaciers with no data from both the `gtd_grids` and `glaciers` vectors and returns the updated vectors.

# Notes
- The function uses `pmap` for parallel processing of glaciers.
- The list of missing glaciers is stored in a JLD2 file located at `params.simulation.working_dir/data/missing_glaciers.jld2`.
- Glaciers with no data are identified and removed based on the condition that all data points in their thickness grid are zero.
"""
function get_glathida!(glaciers::Vector{G}, params::Parameters; force=false) where {G <: Glacier2D}
    gtd_grids = pmap(glacier -> get_glathida_glacier(glacier, params, force), glaciers)

     # Update missing_glaciers list before removing them
    missing_glaciers = load(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
    for (gtd_grid, glacier) in zip(gtd_grids, glaciers)
        if (length(gtd_grid[gtd_grid .!= 0.0]) == 0) && all(glacier.rgi_id .!= missing_glaciers)
            push!(missing_glaciers, glacier.rgi_id)
            @info "Glacier with all data at 0: $(glacier.rgi_id). Updating list of missing glaciers..."
        end
    end
    jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"); missing_glaciers)

   # Apply deletion to both gtd_grids and glaciers using the same set of indices
    indices_to_remove = findall(x -> length(x[x .!= 0.0]) == 0, gtd_grids)
    deleteat!(gtd_grids, indices_to_remove)
    deleteat!(glaciers, indices_to_remove)

    return gtd_grids, glaciers
end

"""
    get_glathida_glacier(glacier::Glacier2D, params::Parameters, force)

Retrieve or generate the glathida glacier grid for a given glacier.

# Arguments
- `glacier::Glacier2D`: The glacier object for which the glathida grid is to be retrieved or generated.
- `params::Parameters`: The parameters object containing simulation settings.
- `force`: A boolean flag indicating whether to force regeneration of the glathida grid even if it already exists.

# Returns
- `gtd_grid`: A 2D array representing the glathida glacier grid.

# Description
This function checks if the glathida glacier grid file (`glathida.h5`) exists in the specified path. If the file exists and `force` is `false`, it reads the grid from the file. Otherwise, it reads the glacier thickness data from a CSV file (`glathida_data.csv`), computes the average thickness for each grid cell, and saves the resulting grid to an HDF5 file (`glathida.h5`).
"""
function get_glathida_glacier(glacier::Glacier2D, params::Parameters, force)
    rgi_path = joinpath(prepro_dir, params.simulation.rgi_paths[glacier.rgi_id])
    gtd_path = joinpath(rgi_path, "glathida.h5")
    if isfile(gtd_path) && !force
        gtd_grid = h5read(gtd_path, "gtd_grid")
    else
        glathida = CSV.File(joinpath(rgi_path, "glathida_data.csv"))
        gtd_grid = zeros(size(glacier.H₀))
        count = zeros(size(glacier.H₀))
        for (thick, i, j) in zip(glathida["thickness"], glathida["i_grid"], glathida["j_grid"])
            count[i,j] += 1
            gtd_grid[i,j] += thick
        end

        gtd_grid .= ifelse.(count .> 0, gtd_grid ./ count, 0.0)

        # Save file
        h5open(joinpath(rgi_path, "glathida.h5"), "w") do file
            write(file, "gtd_grid", gtd_grid)
        end
    end
    return gtd_grid
end

# [End] Glathida Utilities

"""
    filter_missing_glaciers!(rgi_ids::Vector{String}, params::Parameters)

Filter out glaciers that cannot be processed from the given list of RGI IDs.

# Arguments
- `rgi_ids::Vector{String}`: A vector of RGI IDs representing glaciers.
- `params::Parameters`: A `Parameters` object containing simulation parameters.

# Description
This function filters out glaciers from the provided `rgi_ids` list based on two criteria:
1. Glaciers that are marked as level 2 in the RGI statistics CSV file.
2. Glaciers listed in the `missing_glaciers.jld2` file located in the `params.simulation.working_dir` directory.

# Notes
- The RGI statistics CSV file is downloaded from a remote server.
- If the `missing_glaciers.jld2` file is not available, a warning is logged and the function skips this filtering step.
"""
function filter_missing_glaciers!(rgi_ids::Vector{String}, params::Parameters) # TODO: see if this is necessary, otherwise remove

    # Check which glaciers we can actually process
    pathCsv = joinpath(dirname(prepro_dir), "rgi62_stats.csv")
    rgi_stats = CSV.File(pathCsv)

    # Remove level 2 glaciers
    for rgi_id in rgi_ids
        if rgi_stats.Connect[rgi_stats.RGIId .== rgi_id] == 2
            @warn "Filtering glacier $rgi_id..."
            deleteat!(rgi_ids, rgi_ids .== rgi_id)
        end

    end

    try
        missing_glaciers = load(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
        for missing_glacier in missing_glaciers
            deleteat!(rgi_ids, findall(x->x == missing_glacier,rgi_ids))
        end
        @info "Filtering out these glaciers from RGI ID list: $missing_glaciers"
    catch error
        @warn "$error: No missing_glaciers.jld file available. Skipping..."
    end

end

"""
    glacierName(rgi_id::String)
    glacierName(rgi_ids::Vector{String})

Returns the name(s) of one or multiple glaciers based the given RGI ID(s).
It uses the `rgi62_stats.csv` file from OGGM.
"""
function glacierName(rgi_id::String)
    return glacierName([rgi_id])[1]
end
function glacierName(rgi_ids::Vector{String})
    pathCsv = joinpath(dirname(prepro_dir), "rgi62_stats.csv")
    rgi_stats = CSV.File(pathCsv)
    return [glacierName(rgi_id, rgi_stats) for rgi_id in rgi_ids]
end
function glacierName(rgi_id::String, rgi_stats)
    name = rgi_stats.Name[rgi_stats.RGIId .== rgi_id]
    name = replace(name, missing => "")
    if length(name)==0
        @warn "RGI ID $(rgi_id) has no corresponding entry in rgi62_stats."
        name = [""]
    elseif length(name)>1
        @warn "RGI ID $(rgi_id) has multiple corresponding entries in rgi62_stats."
    end
    return name[1]
end

"""
    fillNaN!(A::AbstractArray, fill::Number=zero(eltype(A)))

Replace all `NaN` values in the array `A` with the specified `fill` value.

# Arguments
- `A::AbstractArray`: The array in which `NaN` values will be replaced.
- `fill::Number`: The value to replace `NaN` with. Defaults to `zero(eltype(A))`.
"""
function fillNaN!(A, fill=zero(eltype(A)))
    for i in eachindex(A)
        @inbounds A[i] = ifelse(isnan(A[i]), fill, A[i])
    end
end

"""
    fillNaN(A::AbstractArray, fill::Number=zero(eltype(A)))

Replace all NaN values in the array `A` with the specified `fill` value. 
If no `fill` value is provided, it defaults to the zero value of the element type of `A`.

# Arguments
- `A::AbstractArray`: The input array that may contain NaN values.
- `fill::Number`: The value to replace NaNs with. Defaults to `zero(eltype(A))`.

# Returns
- An array of the same type and shape as `A`, with all NaN values replaced by `fill`.
"""
function fillNaN(A, fill=zero(eltype(A)))
    return @. ifelse(isnan(A), fill, A)
end

"""
    fillZeros!(A::AbstractArray, fill::Number=NaN)

Replace all zero elements in the array `A` with the specified `fill` value.

# Arguments
- `A::AbstractArray`: The array in which to replace zero elements.
- `fill::Number`: The value to replace zero elements with. Defaults to `NaN`.
"""
function fillZeros!(A, fill=NaN)
    for i in eachindex(A)
        @inbounds A[i] = ifelse(iszero(A[i]), fill, A[i])
    end
end

"""
    fillZeros(A::AbstractArray, fill::Number=NaN) -> AbstractArray

Replace all zero elements in the array `A` with the specified `fill` value.

# Arguments
- `A::AbstractArray`: The input array in which zero elements are to be replaced.
- `fill::Number`: The value to replace zero elements with. Defaults to `NaN`.

# Returns
- `AbstractArray`: A new array with zero elements replaced by the `fill` value.
"""
function fillZeros(A, fill=NaN)
    return @. ifelse(iszero(A), fill, A)
end

"""
    parse_proj(proj::String)

Parses the string containing the information of the projection to filter for important information
"+proj=tmerc +lat_0=0 +lon_0=6.985 +k=0.9996 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
"""
function parse_proj(proj::String)
    res = Dict()
    ℓ = split(proj, (' ', '+', '='))
    ℓ = ℓ[ℓ .!= ""]
    for (i, key) in enumerate(ℓ)
        if key ∈ ["lat_0", "lon_0", "k", "x_0", "y_0", "zone"]
            res[key] = parse(Float64, ℓ[i+1])
        end
    end
    return res
end

"""
    UTMercator(x::F, y::F; k=0.9996, cenlon=0.0, cenlat=0.0, x0=0.0, y0=0.0, zone::Union{Nothing, Int}=nothing, hemisphere=nothing) where {F <: AbstractFloat}

Transverse Mercator Projection.
This function reprojects northing/easting coordinates into latitude/longitude.

Keyword arguments
=================
    - `k`: scale factor of the projection
    - `cenlon`: Central longitude used in the projection
    - `cenlat`: Central latitude used in the projection
    - `x0`: Shift in easting
    - `y0`: Shift in northing
    - `zone` : Zone of the projection
    - `hemisphere`: Either :north or :south
"""
function UTMercator(x::F, y::F; k=0.9996, cenlon=0.0, cenlat=0.0, x0=0.0, y0=0.0, zone::Union{Nothing, Int}=nothing, hemisphere=nothing) where {F <: AbstractFloat}

    if !isnothing(zone)
        @assert !isnothing(hemisphere) "When zone is provided, hemisphere should also be defined. It can be either :north or :south"
        projection = CoordRefSystems.utm(hemisphere, zone; datum = WGS84Latest)(x,y)
    else
        # Convert to right units
        lonₒ = cenlon * 1.0°
        latₒ = cenlat * 1.0°
        xₒ = x0 * 1.0m
        yₒ = y0 * 1.0m
        # Define shift in new coordinate system
        S = CoordRefSystems.Shift(; lonₒ, xₒ, yₒ)
        # Define custom projection
        projection = TransverseMercator{k, latₒ, WGS84Latest, S}(x, y)
    end

    return convert(LatLon, projection)
end

"""
    ReverseUTMercator(x::F, y::F; k=0.9996, cenlon=0.0, cenlat=0.0, x0=0.0, y0=0.0, zone::Union{Nothing, Int}=nothing, hemisphere=nothing) where {F <: AbstractFloat}

Transverse Mercator Projection.
This function reprojects latitude/longitude into northing/easting coordinates.

Keyword arguments
=================
    - `k`: scale factor of the projection
    - `cenlon`: Central longitude used in the projection
    - `cenlat`: Central latitude used in the projection
    - `x0`: Shift in easting
    - `y0`: Shift in northing
    - `zone` : Zone of the projection
    - `hemisphere`: Either :north or :south
"""
function ReverseUTMercator(lat::F, lon::F; k=0.9996, cenlon=0.0, cenlat=0.0, x0=0.0, y0=0.0, zone::Union{Nothing, Int}=nothing, hemisphere=nothing) where {F <: AbstractFloat}

    if !isnothing(zone)
        @assert !isnothing(hemisphere) "When zone is provided, hemisphere should also be defined. It can be either :north or :south"
        projection = CoordRefSystems.utm(hemisphere, zone; datum = WGS84Latest)
    else
        # Convert to right units
        lonₒ = cenlon * 1.0°
        latₒ = cenlat * 1.0°
        xₒ = x0 * 1.0m
        yₒ = y0 * 1.0m
        # Define shift in new coordinate system
        S = CoordRefSystems.Shift(; lonₒ, xₒ, yₒ)
        # Define custom projection
        projection = TransverseMercator{k, latₒ, WGS84Latest, S}
    end

    latlon = LatLon(lat, lon)
    return convert(projection, latlon)
end

"""
    smooth!(A)

Smooths the interior of a 2D array `A` using a simple averaging method. The function modifies the array `A` in place.

# Arguments
- `A::AbstractMatrix`: A 2D array to be smoothed.

# Details
The function updates the interior elements of `A` (excluding the boundary elements) by adding a weighted average of the second differences along both dimensions. The boundary elements are then set to the values of their nearest interior neighbors to maintain the boundary conditions.
"""
@views function smooth!(A)
    A[2:end-1,2:end-1] .= A[2:end-1,2:end-1] .+ 1.0./4.1.*(diff(diff(A[:,2:end-1], dims=1), dims=1) .+ diff(diff(A[2:end-1,:], dims=2), dims=2))
    A[1,:]=A[2,:]; A[end,:]=A[end-1,:]; A[:,1]=A[:,2]; A[:,end]=A[:,end-1]
end

# function smooth(A)
#     A_smooth = A[2:end-1,2:end-1] .+ 1.0./4.1.*(diff(diff(A[:,2:end-1], dims=1), dims=1) .+ diff(diff(A[2:end-1,:], dims=2), dims=2))
#     @tullio A_smooth_pad[i,j] := A_smooth[pad(i-1,1,1),pad(j-1,1,1)] # Fill borders 
#     return A_smooth_pad
# end

"""
    is_in_glacier(A::Matrix{F}, distance::I) where {I <: Integer, F <: AbstractFloat}

Return a matrix with booleans indicating if a given pixel is at distance at least
`distance` in the set of non zero values of the matrix. This usually allows
discarding the border pixels of a glacier.
A positive value of `distance`` indicates a measurement from inside the glacier, while a
negative `distance`` indicates one from outside.

Arguments:
- `A::Matrix{F}`: Matrix from which to compute the matrix of booleans.
- `distance::I`: Distance to the border, computed as the number of pixels we need
    to move from within the glacier to find a pixel with value zero.
"""
function is_in_glacier(A::Matrix{F}, distance::I) where {I <: Integer, F <: AbstractFloat}
    B = convert.(F, (A .!= 0))
    # Reverse values in case we want distance from outside the border
    if distance < 0
        distance = -distance
        B .= 1.0 .- B
    end
    for i in 1:distance
        # We cannot use in-place affectation because this function is differentiated by Zygote in ODINN
        B = min.(
            B,
            circshift(B, (1,0)),
            circshift(B, (-1,0)),
            circshift(B, (0,1)),
            circshift(B, (0,-1))
            )
    end
    B_bool = B .> 0.001
    if distance >= 0
        return B_bool
    else
        return .!B_bool
    end
end

"""
    block_average_pad_edge(mat::Matrix{F}, n::Int) where {F <: AbstractFloat}

Downsamples a matrix by averaging `n x n` blocks, using edge-replication padding
when the matrix dimensions are not divisible by `n`.
Edge padding replicates the last row/column values to expand the matrix so that both
dimensions are divisible by `n`.
Returns a matrix of averaged values with size `(ceil(Int, X/n), ceil(Int, Y/n))`.

Arguments
- `mat::Matrix{F}`: Input 2D matrix.
- `n::Int`: Block size for downsampling.
"""
function block_average_pad_edge(mat::Matrix{F}, n::Int) where {F <: AbstractFloat}
    X, Y = size(mat)
    new_X = ceil(Int, X / n) * n
    new_Y = ceil(Int, Y / n) * n

    # Create padded matrix filled with edge values
    padded = similar(mat, new_X, new_Y)

    # Fill original data
    padded[1:X, 1:Y] .= mat

    # Pad bottom rows with last row
    if new_X > X
        for i in X+1:new_X
            padded[i, 1:Y] .= mat[end, :]
        end
    end

    # Pad right columns with last column
    if new_Y > Y
        for j in Y+1:new_Y
            padded[1:X, j] .= mat[:, end]
        end
    end

    # Fill bottom-right corner (if both X and Y were padded)
    if new_X > X && new_Y > Y
        for i in X+1:new_X
            for j in Y+1:new_Y
                padded[i, j] = mat[end, end]
            end
        end
    end

    return block_average(padded, n)
end

"""
    block_average(mat::Matrix{F}, n::Int) where {F <: AbstractFloat}

Downsamples a matrix by averaging non-overlapping `n x n` blocks.
Returns a matrix of the block-averaged values with size `(div(X, n), div(Y, n))`
where `(X, Y) = size(mat)`.

Arguments
- `mat::Matrix{F}`: Input 2D matrix.
- `n::Int`: Block size for downsampling. Both matrix dimensions must be divisible by `n`.
"""
function block_average(mat::Matrix{F}, n::Int) where {F <: AbstractFloat}
    X, Y = size(mat)
    @assert X % n == 0 && Y % n == 0 "Matrix dimensions are $(size(mat)) but they are not divisible by n=$n"

    A, B = div(X, n), div(Y, n)
    reshaped = reshape(mat, n, A, n, B)
    permuted = permutedims(reshaped, (2, 4, 1, 3))  # (A, B, n, n)
    mean_blocks = mean(permuted, dims=(3, 4))
    return dropdims(mean_blocks, dims=(3, 4))
end
