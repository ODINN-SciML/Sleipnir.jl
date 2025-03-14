
export initialize_glaciers

###############################################
############  FUNCTIONS   #####################
###############################################

"""
    initialize_glaciers(rgi_ids::Vector{String}, params::Parameters; test=false)

Initialize glaciers based on provided RGI IDs and parameters.

# Arguments
- `rgi_ids::Vector{String}`: A vector of RGI IDs representing the glaciers to be initialized.
- `params::Parameters`: A `Parameters` object containing simulation parameters.
- `test::Bool`: An optional boolean flag indicating whether to run in test mode. Default is `false`.

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
function initialize_glaciers(rgi_ids::Vector{String}, params::Parameters; test=false)

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
        map((rgi_id) -> generate_raw_climate_files(rgi_id, params.simulation), rgi_ids) # avoid GitHub CI issue
        # glaciers::Vector{Glacier2D} = map((rgi_id) -> initialize_glacier(rgi_id, params; smoothing=false, test=test), rgi_ids)
    else
        pmap((rgi_id) -> generate_raw_climate_files(rgi_id, params.simulation), rgi_ids)
        # glaciers::Vector{Glacier2D} = pmap((rgi_id) -> initialize_glacier(rgi_id, params; smoothing=false, test=test), rgi_ids)
    end
        
    glaciers::Vector{Glacier2D} = pmap((rgi_id) -> initialize_glacier(rgi_id, params; smoothing=false, test=test), rgi_ids)

    if params.simulation.use_glathida_data == true
        
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
    initialize_glacier(rgi_id::String, parameters::Parameters; smoothing=false, test=false)

Initialize a glacier with the given RGI ID and parameters.

# Arguments
- `rgi_id::String`: The RGI (Randolph Glacier Inventory) ID of the glacier.
- `parameters::Parameters`: A struct containing various parameters required for initialization.
- `smoothing::Bool`: Optional. If `true`, apply smoothing to the initial topography. Default is `false`.
- `test::Bool`: Optional. If `true`, run in test mode. Default is `false`.

# Returns
- `glacier`: An initialized glacier object containing the initial topography and climate data.
"""
function initialize_glacier(rgi_id::String, parameters::Parameters; smoothing=false, test=false)
    # Initialize glacier initial topography
    glacier = initialize_glacier_data(rgi_id, parameters; smoothing=smoothing, test=test)

    # Initialize glacier climate
    initialize_glacier_climate!(glacier, parameters)

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
    initialize_glacier_data(rgi_id::String, params::Parameters; smoothing=false, test=false)

Initialize glacier data for a given RGI ID and parameters.

# Arguments
- `rgi_id::String`: The RGI ID of the glacier.
- `params::Parameters`: A `Parameters` object containing simulation parameters.
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
function initialize_glacier_data(rgi_id::String, params::Parameters; smoothing=false, test=false)
    # Load glacier gridded data
    F = Sleipnir.Float
    rgi_path = joinpath(prepro_dir, params.simulation.rgi_paths[rgi_id])
    glacier_gd = RasterStack(joinpath(rgi_path, "gridded_data.nc"))
    if Sleipnir.doublePrec
        glacier_gd = convertRasterStackToFloat64(glacier_gd)
    end
    glacier_grid = JSON.parsefile(joinpath(rgi_path, "glacier_grid.json"))
    # println("Using $ice_thickness_source for initial state")
    # Retrieve initial conditions from OGGM
    # initial ice thickness conditions for forward model
    if params.simulation.ice_thickness_source == "Millan22" && params.simulation.velocities
        H₀ = F.(reverse((ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_ice_thickness.data, 0.0)), dims=2)) # all matrices are reversed
    elseif params.simulation.ice_thickness_source == "Farinotti19"
        H₀ = F.(reverse((ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.consensus_ice_thickness.data, 0.0)), dims=2)) # all matrices are reversed
    end
    fillNaN!(H₀) # Fill NaNs with 0s to have real boundary conditions
    if smoothing
        println("Smoothing is being applied to initial condition.")
        smooth!(H₀)  # Smooth initial ice thickness to help the solver
    end

    try
        # We filter glacier borders in high elevations to avoid overflow problems
        dist_border::Matrix{Sleipnir.Float} = reverse(glacier_gd.dis_from_border.data, dims=2) # matrix needs to be reversed
        
            # H_mask = (dist_border .< 20.0) .&& (S .> maximum(S)*0.7)
            # H₀[H_mask] .= 0.0

        # Mercator Projection 
        params_projection = parse_proj(glacier_grid["proj"])
        transform(X,Y) = UTMercator(X, Y; k=params_projection["k"], cenlon=params_projection["lon_0"], cenlat=params_projection["lat_0"], 
                        x0=params_projection["x_0"], y0=params_projection["y_0"])
        easting = dims(glacier_gd, 1).val
        northing = dims(glacier_gd, 2).val
        latitudes = map(x -> x.lat.val, transform.(Ref(mean(easting)), northing))
        longitudes = map(x -> x.lon.val, transform.(easting, Ref(mean(northing))))
        x0y0 = transform(glacier_grid["x0y0"][1], glacier_grid["x0y0"][2])
        cenlon::Sleipnir.Float = x0y0.lon.val
        cenlat::Sleipnir.Float = x0y0.lat.val
        if maximum(abs.(latitudes)) > 80
            @warn "Mercator projection can fail in high-latitude regions. You glacier includes latitudes larger than 80°."
        end

        B = reverse(glacier_gd.topo.data, dims=2) .- H₀ # bedrock (matrix also needs to be reversed)
        
        Coords = Dict{String,Vector{Float64}}("lon"=> longitudes, "lat"=> latitudes)
        S::Matrix{Sleipnir.Float} = reverse(glacier_gd.topo.data, dims=2)
        #smooth!(S)

        if params.simulation.velocities
            # All matrices need to be reversed
            V::Matrix{Sleipnir.Float} = reverse(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_v.data, 0.0), dims=2)
            Vx::Matrix{Sleipnir.Float} = reverse(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_vx.data, 0.0), dims=2)
            Vy::Matrix{Sleipnir.Float} = reverse(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_vy.data, 0.0), dims=2)
            fillNaN!(V)
            fillNaN!(Vx)
            fillNaN!(Vy)
        else
            V = zeros(F, size(H₀))
            Vx = zeros(F, size(H₀))
            Vy = zeros(F, size(H₀))
        end
        nx = glacier_grid["nxny"][1]
        ny = glacier_grid["nxny"][2]
        Δx::Sleipnir.Float = abs.(glacier_grid["dxdy"][1])
        Δy::Sleipnir.Float = abs.(glacier_grid["dxdy"][2])
        slope::Matrix{Sleipnir.Float} = glacier_gd.slope.data

        # We initialize the Glacier with all the initial topographical 
        glacier = Glacier2D(rgi_id = rgi_id, 
                          climate=nothing, 
                          H₀ = H₀, S = S, B = B, V = V, Vx = Vx, Vy = Vy,
                          A = Sleipnir.Float(4e-17), C = Sleipnir.Float(0.0), n = Sleipnir.Float(3.0),
                          slope = slope, dist_border = dist_border,
                          Coords = Coords, Δx=Δx, Δy=Δy, nx=nx, ny=ny,
                          cenlon = cenlon, cenlat = cenlat)
        return glacier

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
    get_glathida!(glaciers::Vector{Glacier2D}, params::Parameters; force=false)

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
function get_glathida!(glaciers::Vector{Glacier2D}, params::Parameters; force=false)
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
    filter_missing_glaciers!(glaciers::Vector{Glacier2D}, params::Parameters)

Filters out glaciers from the provided `glaciers` vector that are marked as missing in the task log or in a previously saved file.

# Arguments
- `glaciers::Vector{Glacier2D}`: A vector of `Glacier2D` objects to be filtered.
- `params::Parameters`: A `Parameters` object containing simulation parameters.

# Returns
- `missing_glaciers::Vector{String}`: A vector of glacier IDs that were filtered out.

# Details
The function reads a task log CSV file from the working directory specified in `params`. It then determines which glaciers are missing based on the task log and additional conditions specified in `params`. If a previously saved file of missing glaciers exists, it loads and merges the missing glaciers from that file. Finally, it removes the missing glaciers from the `glaciers` vector and saves the updated list of missing glaciers to a file.
"""
function filter_missing_glaciers!(glaciers::Vector{Glacier2D}, params::Parameters)
    task_log = CSV.File(joinpath(params.simulation.working_dir, "task_log.csv"))
    if params.simulation.velocities & params.simulation.use_glathida_data
        glacier_filter = (task_log.velocity_to_gdir .!= "SUCCESS") .&& (task_log.gridded_attributes .!= "SUCCESS") .&& (task_log.thickness_to_gdir .!= "SUCCESS")
    elseif params.simulation.use_glathida_data
        glacier_filter = (task_log.gridded_attributes .!= "SUCCESS") .&& (task_log.thickness_to_gdir .!= "SUCCESS")
    else
        glacier_filter = (task_log.gridded_attributes .!= "SUCCESS")
    end
    
    glacier_ids = Vector{String}([])

    for id in task_log["index"]
        push!(glacier_ids, id)
    end
    missing_glaciers = glacier_ids[glacier_filter]

    try
        missing_glaciers_old = load(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
        for missing_glacier in missing_glaciers_old
            @show missing_glacier
            if all(missing_glacier .!= missing_glaciers) # if the glacier is already not present, let's add it
                push!(missing_glaciers, missing_glacier)
            end
        end
    catch error
        @warn "$error: No missing_glaciers.jld file available. Skipping..."
    end

    for id in missing_glaciers
        deleteat!(glaciers, findall(x->x.rgi_id==id, glaciers))
    end
    
    # Save missing glaciers in a file
    jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"); missing_glaciers)
    #@warn "Filtering out these glaciers from gdir list: $missing_glaciers"
    
    return missing_glaciers
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
        if key ∈ ["lat_0", "lon_0", "k", "x_0", "y_0"]
            res[key] = parse(Float64, ℓ[i+1])
        end
    end
    return res
end

"""
    UTMercator(x::F, y::F; k=0.9996, cenlon=0.0, cenlat=0.0, x0=0.0, y0=0.0)

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
function UTMercator(x::F, y::F; k=0.9996, cenlon=0.0, cenlat=0.0, x0=0.0, y0=0.0, zone::Union{Nothing, Int}=nothing, hemisphere=:north) where {F <: AbstractFloat}
  
    if !isnothing(zone)
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

