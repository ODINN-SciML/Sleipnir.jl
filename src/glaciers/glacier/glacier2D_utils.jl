
export initialize_glaciers

###############################################
############  FUNCTIONS   #####################
###############################################

"""
    initialize_glaciers(rgi_ids::Vector{String}, params::Parameters; test=false)

Initialize multiple `Glacier`s based on a list of RGI IDs and on parameters.

Keyword arguments
=================
    - `rgi_ids`: List of RGI IDs of glaciers
    - `params`: `Parameters` object to be passed
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
    initialize_glacier(rgi_id::String, parameters::Parameters; smoothing=false, velocities=true)

Initialize a single `Glacier`s, including its `Climate`, based on a `rgi_id` and timestepping arguments.

Keyword arguments
=================
    - `rgi_id`: Glacier RGI ID
    - `parameters`: Parameters including the physical and simulation ones
    - `smoothing` Flag determining if smoothing needs to be applied to the surface elevation and ice thickness.
    - `velocities` Flag determining if the ice surface velocities need to be retrieved.
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
    initialize_glacier(rgi_id::String, params::Parameters; smoothing=false, velocities=true)

Retrieves the initial glacier geometry (bedrock + ice thickness) for a glacier with other necessary data (e.g. grid size and ice surface velocities).
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


function filter_missing_glaciers!(rgi_ids::Vector{String}, params::Parameters) # TODO: see if this is necessary, otherwise remove

    # Check which glaciers we can actually process
    pathCsv = Downloads.download("https://cluster.klima.uni-bremen.de/~oggm/rgi/rgi62_stats.csv")
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
fillNaN!(x, fill)

Convert empty matrix grid cells into fill value
"""
function fillNaN!(A, fill=zero(eltype(A)))
    for i in eachindex(A)
        @inbounds A[i] = ifelse(isnan(A[i]), fill, A[i])
    end
end

function fillNaN(A, fill=zero(eltype(A)))
    return @. ifelse(isnan(A), fill, A)
end

function fillZeros!(A, fill=NaN)
    for i in eachindex(A)
        @inbounds A[i] = ifelse(iszero(A[i]), fill, A[i])
    end
end

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

Smooth data contained in a matrix with one time step (CFL) of diffusion.
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

