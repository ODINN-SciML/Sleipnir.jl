
export initialize_glaciers
export glacierName

###############################################################
########  GLACIER INITIALIZATION  ############################
###############################################################

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
        velocityDatacubes::Union{
            Dict{String, String},
            Dict{String, Vector{String}},
            Dict{String, <:RasterStack},
            Dict{String, Vector{<:RasterStack}}
        } = Dict{String, String}(),
        velocityFlag::Union{String, <: RasterStack, Nothing} = nothing
)
    if params.simulation.catch_errors
        # Generate missing glaciers file
        missing_glaciers_path = joinpath(params.simulation.working_dir, "data")
        if !isdir(missing_glaciers_path)
            mkdir(missing_glaciers_path)
        end

        if !isfile(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))
            missing_glaciers = Vector([])
            jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2");
                missing_glaciers)
        end
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
        (rgi_id) -> initialize_glacier(
            rgi_id,
            params;
            velocityDatacubes = velocityDatacubes,
            velocityFlag = velocityFlag
        ),
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
        velocityDatacubes::Union{
            Dict{String, String},
            Dict{String, Vector{String}},
            Dict{String, <: RasterStack},
            Dict{String, Vector{<:RasterStack}}
        } = Dict{String, String}(),
        velocityFlag::Union{String, <: RasterStack, Nothing} = nothing
)
    # Build glacier and its associated climate
    glacier = Glacier2D(rgi_id, parameters; masking = masking, smoothing = smoothing)

    if get(velocityDatacubes, glacier.rgi_id, "") != ""
        mapping = parameters.simulation.mapping
        datacubes = velocityDatacubes[glacier.rgi_id]
        refVelocity = if typeof(datacubes) <: Vector
            refVelocities = map(
                datacube -> initialize_surfacevelocitydata(
                    datacube;
                    glacier = glacier,
                    mapping = mapping,
                    flag = velocityFlag
                ),
                datacubes
            )
            # Combine velocity data with unique datetime
            combine_velocity_data(refVelocities; merge = true)
        else
            # Note: we don't merge data here based on unique date
            initialize_surfacevelocitydata(
                datacubes;
                glacier = glacier,
                mapping = mapping,
                flag = velocityFlag
            )
        end
        # Rebuild glacier since we cannot change type of `glacier.velocityData`
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
        metadata = metadata(rs)
    )
end

function _build_glacier(params, glacier_gd, masking, masking_loss, glacier_grid, H₀, rgi_id)
    # We filter glacier borders in high elevations to avoid overflow problems
    dist_border::Matrix{Sleipnir.Float} = glacier_gd.dis_from_border.data
    if params.simulation.gridScalingFactor > 1
        # Note: this is not mathematically correct and this should be fixed in the future, however since this option is used only in the tests it isn't critical
        dist_border = block_average_pad_edge(
            dist_border, params.simulation.gridScalingFactor)
    end

    # Define mask where ice is constrained to be zero
    mask = @match masking begin
        ::Nothing => falses(size(H₀)...)
        ::Int => is_in_glacier(H₀, -masking)
        ::BitMatrix => masking
    end

    # Define mask where losses used for inversion will be evaluated
    mask_loss = @match masking_loss begin
        ::Nothing => falses(size(H₀)...)
        ::Int => is_in_glacier(H₀, masking_loss)
        ::BitMatrix => masking_loss
    end

    nx, ny = params.simulation.gridScalingFactor > 1 ? size(H₀) : glacier_grid["nxny"]

    # Mercator Projection
    params_projection::Dict{String, Float64} = parse_proj(glacier_grid["proj"])
    transform(X,
        Y) = UTMercator(
        X, Y;
        k = params_projection["k"],
        cenlon = params_projection["lon_0"], cenlat = params_projection["lat_0"],
        x0 = params_projection["x_0"], y0 = params_projection["y_0"]
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

    Coords = Dict{String, Vector{Float64}}("lon" => longitudes, "lat" => latitudes)

    if params.simulation.use_velocities
        V = ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_v.data, 0.0)
        Vx = ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_vx.data, 0.0)
        Vy = ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_vy.data, 0.0)
        fillNaN!(V)
        fillNaN!(Vx)
        fillNaN!(Vy)
    else
        V = zeros(Sleipnir.Float, size(H₀))
        Vx = zeros(Sleipnir.Float, size(H₀))
        Vy = zeros(Sleipnir.Float, size(H₀))
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
        A = Sleipnir.Float(4e-17), C = Sleipnir.Float(0.0),
        n = Sleipnir.Float(3.0), p = Sleipnir.Float(3.0), q = Sleipnir.Float(2.0),
        slope = slope, dist_border = dist_border,
        mask = mask, mask_loss = mask_loss,
        Coords = Coords, Δx = Δx, Δy = Δy, nx = nx, ny = ny,
        cenlon = cenlon, cenlat = cenlat,
        params_projection = params_projection,
        dhdtData = _default_hugonnet_dhdt(rgi_id),
        geodetic_MB = begin
            dhdt_data = _default_hugonnet_dhdt(rgi_id)
            isnothing(dhdt_data) ? Sleipnir.Float(NaN) : Sleipnir.Float(dhdt_data.dhdt)
        end,
        geodetic_MB_uncertainty = _default_hugonnet_mb_uncertainty(rgi_id)
    )
end

"""
    Glacier2D(
        rgi_id::String,
        params::Parameters;
        masking::Union{Int, Nothing, BitMatrix} = 2,
        smoothing=false
    )

Build glacier object for a given RGI ID and parameters.

# Arguments

  - `rgi_id::String`: The RGI ID of the glacier.

  - `params::Parameters`: A `Parameters` object containing simulation parameters.

  - `masking::Union{Int, Nothing, BitMatrix}`: Type of mask applied to the glacier to determine regions with no ice.

      + When `masking` is an `Int`, the mask is based on the initial ice thickness `H₀` and it is set to true for
        pixels outside at a distance of the glacier borders greater than the value of `masking`.
      + When `masking` is set to `nothing`, the mask is set to a `BitMatrix` full of falses.
      + When `masking` is a `BitMatrix`, this matrix is used for the mask.
        Defaults to `2`.

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
        masking_loss::Union{Int, Nothing, BitMatrix} = 0,
        smoothing = false
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
    if params.simulation.ice_thickness_source == :Millan22 &&
       params.simulation.use_velocities
        H₀ = F.(ifelse.(
            glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_ice_thickness.data, 0.0))
    elseif params.simulation.ice_thickness_source == :Farinotti19
        H₀ = F.(ifelse.(glacier_gd.glacier_mask.data .== 1,
            glacier_gd.consensus_ice_thickness.data, 0.0))
    end
    fillNaN!(H₀) # Fill NaNs with 0s to have real boundary conditions
    if smoothing
        println("Smoothing is being applied to initial condition.")
        smooth!(H₀)  # Smooth initial ice thickness to help the solver
    end
    if params.simulation.gridScalingFactor > 1
        H₀ = block_average_pad_edge(H₀, params.simulation.gridScalingFactor)
    end

    if params.simulation.catch_errors
        try
            return _build_glacier(
                params, glacier_gd, masking, masking_loss, glacier_grid, H₀, rgi_id)

        catch error
            @show error
            missing_glaciers = load(joinpath(
                params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
            push!(missing_glaciers, rgi_id)
            jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2");
                missing_glaciers)
            @warn "Glacier without data: $rgi_id. Updating list of missing glaciers. Please try again."
        end
    else
        return _build_glacier(
            params, glacier_gd, masking, masking_loss, glacier_grid, H₀, rgi_id)
    end
end

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

# TODO: see if this is necessary, otherwise remove

  - The RGI statistics CSV file is downloaded from a remote server.
  - If the `missing_glaciers.jld2` file is not available, a warning is logged and the function skips this filtering step.    # Check which glaciers we can actually process # TODO: see if this is necessary, otherwise remove
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

    if params.simulation.catch_errors
        try
            missing_glaciers = load(joinpath(
                params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
            for missing_glacier in missing_glaciers
                deleteat!(rgi_ids, findall(x->x == missing_glacier, rgi_ids))
            end
            if length(missing_glaciers) > 0
                @info "Filtering out these glaciers from RGI ID list: $missing_glaciers"
            end
        catch error
            @warn "$error: No missing_glaciers.jld file available. Skipping..."
        end
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
function fillNaN!(A, fill = zero(eltype(A)))
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
function fillNaN(A, fill = zero(eltype(A)))
    return @. ifelse(isnan(A), fill, A)
end

"""
    fillZeros!(A::AbstractArray, fill::Number=NaN)

Replace all zero elements in the array `A` with the specified `fill` value.

# Arguments

  - `A::AbstractArray`: The array in which to replace zero elements.
  - `fill::Number`: The value to replace zero elements with. Defaults to `NaN`.
"""
function fillZeros!(A, fill = NaN)
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
function fillZeros(A, fill = NaN)
    return @. ifelse(iszero(A), fill, A)
end

