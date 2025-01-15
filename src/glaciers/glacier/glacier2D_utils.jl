
export initialize_glaciers

###############################################
############  FUNCTIONS   #####################
###############################################

"""
    initialize_glaciers(rgi_ids::Vector{String}, params::Parameters; velocities=true)

Initialize multiple `Glacier`s based on a list of RGI IDs, a º span for a simulation and step.
    
Keyword arguments
=================
    - `rgi_ids`: List of RGI IDs of glaciers
    - `tspan`: Tuple specifying the initial and final year of the simulation
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
    else
        pmap((rgi_id) -> generate_raw_climate_files(rgi_id, params.simulation), rgi_ids)
    end
    
    glaciers = pmap((rgi_id) -> initialize_glacier(rgi_id, params; smoothing=false, test=test), rgi_ids)
    
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
    initialize_glacier(gdir::Py, tspan, step; smoothing=false, velocities=true)

Initialize a single `Glacier`s, including its `Climate`, based on a `gdir` and timestepping arguments.
    
Keyword arguments
=================
    - `gdir`: Glacier directory
    - `tspan`: Tuple specifying the initial and final year of the simulation
    - `step`: Step in years for the surface mass balance processing
    - `smoothing` Flag determining if smoothing needs to be applied to the surface elevation and ice thickness.
    - `velocities` Flag determining if the ice surface velocities need to be retrieved.
"""
function initialize_glacier(rgi_id::String, parameters::Parameters; smoothing=false, test=false)
    # Initialize glacier initial topography
    glacier = initialize_glacier_data(rgi_id, parameters; smoothing=smoothing, test=test)

    # Initialize glacier climate
    initialize_glacier_climate!(glacier, parameters)

    if test
        glacier.gdir = nothing
        glacier.S_coords = nothing
    end

    return glacier
end

"""
    initialize_glacier(gdir::Py; smoothing=false, velocities=true)

Retrieves the initial glacier geometry (bedrock + ice thickness) for a glacier with other necessary data (e.g. grid size and ice surface velocities).
"""
function initialize_glacier_data(rgi_id::String, params::Parameters; smoothing=false, test=false)
    # Load glacier gridded data
    F = params.simulation.float_type
    I = params.simulation.int_type
    rgi_path = params.simulation.rgi_path[rgi_id]
    glacier_gd = RasterStack(joinpath(rgi_path, "gridded_data.nc"))
    glacier_grid = JSON.parsefile(joinpath(rgi_path, "glacier_grid.json"))
    # println("Using $ice_thickness_source for initial state")
    # Retrieve initial conditions from OGGM
    # initial ice thickness conditions for forward model
    if params.OGGM.ice_thickness_source == "Millan22" && params.simulation.velocities
        H₀ = F.(ifelse.(Matrix,glacier_gd.glacier_mask.data .== 1, Matrix,glacier_gd.millan_ice_thickness.data, 0.0))
    elseif params.OGGM.ice_thickness_source == "Farinotti19"
        H₀ = F.(ifelse.(Matrix,glacier_gd.glacier_mask.data .== 1, Matrix,glacier_gd.consensus_ice_thickness.data, 0.0))
    end
    fillNaN!(H₀) # Fill NaNs with 0s to have real boundary conditions
    if smoothing 
        println("Smoothing is being applied to initial condition.")
        smooth!(H₀)  # Smooth initial ice thickness to help the solver
    end

    # # Create path for simulation results
    # gdir_path = dirname(pyconvert(String, gdir.get_filepath("dem")))
    # if !isdir(gdir_path)
    #     mkdir(gdir_path)
    # end

    try
        # We filter glacier borders in high elevations to avoid overflow problems
        dist_border = glacier_gd.dis_from_border.data
        
            # H_mask = (dist_border .< 20.0) .&& (S .> maximum(S)*0.7)
            # H₀[H_mask] .= 0.0

        B = glacier_gd.topo.data .- H₀ # bedrock
        
        S_coords = Dict{"x"=> glacier_gd.topo.dims[1], "y"=> glacier_gd.topo.dims[2]}
        S = glacier_gd.topo.data
        #smooth!(S)
        
        if params.simulation.velocities
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
        nx = glacier_grid["nxny"][1] 
        ny = glacier_grid["nxny"][2] 
        Δx = abs.(glacier_grid["dxdy"][1])
        Δy = abs.(glacier_grid["dxdy"][2])
        slope = glacier_gd.slope.data

        # We initialize the Glacier with all the initial topographical 
        glacier = Glacier2D(rgi_id = rgi_id, 
                          climate=nothing, 
                          H₀ = H₀, S = S, B = B, V = V, Vx = Vx, Vy = Vy,
                          A = 4e-17, C = 0.0, n = 3.0,
                          slope = slope, dist_border = dist_border,
                          S_coords = S_coords, Δx=Δx, Δy=Δy, nx=nx, ny=ny,
                          cenlon = glacier_grid["x0y0"][1] , cenlat = glacier_grid["x0y0"][2])
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

   # Apply deletion to both gtd_grids and gdirs using the same set of indices
    indices_to_remove = findall(x -> length(x[x .!= 0.0]) == 0, gtd_grids)
    deleteat!(gtd_grids, indices_to_remove)
    deleteat!(glaciers, indices_to_remove)
   
    return gtd_grids, glaciers
end

function get_glathida_glacier(glacier::Glacier2D, params::Parameters, force)
    gtd_path = joinpath(gdir.dir, "glathida.h5")
    if isfile(gtd_path) && !force
        gtd_grid = h5read(gtd_path, "gtd_grid")
    else
        glathida = CSV.File(joinpath(params.simulation.rgi_path, "glathida.csv")) 
        gtd_grid = zeros(size(glacier.H₀))
        count = zeros(size(glacier.H₀))
        for (thick, i, j) in zip(glathida["elevation"], glathida["i_grid"], glathida["j_grid"])
            count[i,j] += 1
            gtd_grid[i,j] += thick
        end

        gtd_grid .= ifelse.(count > 0, gtd_grid ./ count, 0.0)
        
        # Save file 
        h5open(joinpath(gdir.dir, "glathida.h5"), "w") do file
            write(file, "gtd_grid", gtd_grid)  
        end
    end
    return gtd_grid
end

function get_glathida_path_and_IDs()
    gtd_file = Downloads.download("https://cluster.klima.uni-bremen.de/~oggm/glathida/glathida-v3.1.0/data/TTT_per_rgi_id.h5")
    glathida = h5open(gtd_file, "r")
    rgi_ids = keys(glathida)
    rgi_ids = String[id[2:end] for id in rgi_ids]
    return gtd_file, rgi_ids
end
# [End] Glathida Utilities


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

function filter_missing_glaciers!(rgi_ids::Vector{String}, params::Parameters) # TODO: see if this is necessary and convert to Julia

    # Check which glaciers we can actually process
    rgi_stats = pd[].read_csv(utils[].file_downloader("https://cluster.klima.uni-bremen.de/~oggm/rgi/rgi62_stats.csv"), index_col=0)
    # rgi_stats = rgi_stats.loc[rgi_ids]

    # if any(rgi_stats.Connect .== 2)
    #     @warn "You have some level 2 glaciers... Removing..."
    #     rgi_ids = [rgi_stats.loc[rgi_stats.Connect .!= 2].index]
    # end

    indices = [rgi_stats.index...]
    for rgi_id in rgi_ids
        if PyList(rgi_stats.Connect.values[indices .== rgi_id]) == 2
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

