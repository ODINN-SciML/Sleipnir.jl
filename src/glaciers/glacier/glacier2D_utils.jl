
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

    # Initialize glacier directories
    gdirs::Vector{PyObject} = init_gdirs(rgi_ids, params; velocities=params.simulation.velocities)
    
    # Generate raw climate data if necessary
    if params.simulation.test_mode
        map((gdir) -> generate_raw_climate_files(gdir, params.simulation.tspan), gdirs) # avoid GitHub CI issue
    else
        pmap((gdir) -> generate_raw_climate_files(gdir, params.simulation.tspan), gdirs)
    end
    
    glaciers::Vector{Glacier2D} = pmap((gdir) -> initialize_glacier(gdir, params; smoothing=false, test=test), gdirs)
    
    if params.simulation.use_glathida_data == true
        data_glathida, glathida_rgi_ids = get_glathida_path_and_IDs()
        
        # Obtain H_glathida values for the valid RGI IDs
        H_glathida_values, valid_gdirs = get_glathida!(data_glathida, gdirs, params)
        valid_rgi_ids = [gdir.rgi_id for gdir in valid_gdirs]
        
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
    initialize_glacier(gdir::PyObject, tspan, step; smoothing=false, velocities=true)

Initialize a single `Glacier`s, including its `Climate`, based on a `gdir` and timestepping arguments.
    
Keyword arguments
=================
    - `gdir`: Glacier directory
    - `tspan`: Tuple specifying the initial and final year of the simulation
    - `step`: Step in years for the surface mass balance processing
    - `smoothing` Flag determining if smoothing needs to be applied to the surface elevation and ice thickness.
    - `velocities` Flag determining if the ice surface velocities need to be retrieved.
"""
function initialize_glacier(gdir::PyObject, parameters; smoothing=false, test=false)
    # Initialize glacier initial topography
    glacier = initialize_glacier_data(gdir, parameters; smoothing=smoothing, test=test)

    # Initialize glacier climate
    initialize_glacier_climate!(glacier, parameters)

    if test
        glacier.gdir = nothing
        glacier.S_coords = nothing
    end

    return glacier
end

"""
    initialize_glacier(gdir::PyObject; smoothing=false, velocities=true)

Retrieves the initial glacier geometry (bedrock + ice thickness) for a glacier with other necessary data (e.g. grid size and ice surface velocities).
"""
function initialize_glacier_data(gdir::PyObject, params::Parameters; smoothing=false, test=false)
    # Load glacier gridded data
    F = params.simulation.float_type
    I = params.simulation.int_type
    glacier_gd = xr.open_dataset(gdir.get_filepath("gridded_data"))
    # println("Using $ice_thickness_source for initial state")
    # Retrieve initial conditions from OGGM
    # initial ice thickness conditions for forward model
    if params.OGGM.ice_thickness_source == "Millan22" && params.simulation.velocities
        H₀ = F.(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_ice_thickness.data, 0.0))
    elseif params.OGGM.ice_thickness_source == "Farinotti19"
        H₀ = F.(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.consensus_ice_thickness.data, 0.0))
    end
    fillNaN!(H₀) # Fill NaNs with 0s to have real boundary conditions
    if smoothing 
        println("Smoothing is being applied to initial condition.")
        smooth!(H₀)  # Smooth initial ice thickness to help the solver
    end

    # Create path for simulation results
    gdir_path = dirname(gdir.get_filepath("dem"))
    if !isdir(gdir_path)
        mkdir(gdir_path)
    end

    try
        # We filter glacier borders in high elevations to avoid overflow problems
        dist_border::Matrix{F} = F.(glacier_gd.dis_from_border.data)
        
            # H_mask = (dist_border .< 20.0) .&& (S .> maximum(S)*0.7)
            # H₀[H_mask] .= 0.0

        B::Matrix{F} = F.(glacier_gd.topo.data) .- H₀ # bedrock
        S_coords::PyObject = rioxarray.open_rasterio(gdir.get_filepath("dem"))
        #S::Matrix{F} = F.(glacier_gd.topo.data)
        S::Matrix{F} = F.(S_coords[1].values) # surface elevation
        #smooth!(S)
        
        if params.simulation.velocities
            V::Matrix{F} = F.(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_v.data, 0.0))
            Vx::Matrix{F} = F.(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_vx.data, 0.0))
            Vy::Matrix{F} = F.(ifelse.(glacier_gd.glacier_mask.data .== 1, glacier_gd.millan_vy.data, 0.0))
            fillNaN!(V)
            fillNaN!(Vx)
            fillNaN!(Vy)
        else
            V = zeros(F, size(H₀))
            Vx = zeros(F, size(H₀))
            Vy = zeros(F, size(H₀))
        end
        nx = I(glacier_gd.y.size) # glacier extent
        ny = I(glacier_gd.x.size) # really weird, but this is inversed 
        Δx = abs(F(gdir.grid.dx))
        Δy = abs(F.(gdir.grid.dy))
        slope = F.(glacier_gd.slope.data)

        glacier_gd.close() # Release any resources linked to this object
        # We initialize the Glacier with all the initial topographical conditions
        glacier = Glacier2D(rgi_id = gdir.rgi_id, gdir = gdir,
                          climate=nothing, 
                          H₀ = H₀, S = S, B = B, V = V, Vx = Vx, Vy = Vy,
                          A = 4e-17, C = 0.0, n = 3.0,
                          slope = slope, dist_border = dist_border,
                          S_coords = S_coords, Δx=Δx, Δy=Δy, nx=nx, ny=ny,
                          cenlon = gdir.cenlon, cenlat = gdir.cenlat)
        return glacier

    catch error
        @show error  
        missing_glaciers = load(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
        push!(missing_glaciers, gdir.rgi_id)
        jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"); missing_glaciers)
        glacier_gd.close() # Release any resources linked to this object
        @warn "Glacier without data: $(gdir.rgi_id). Updating list of missing glaciers. Please try again."
    end
end

"""
    init_gdirs(rgi_ids; force=false)

Initializes Glacier Directories using OGGM. Wrapper function calling `init_gdirs_scratch(rgi_ids)`.
"""
function init_gdirs(rgi_ids::Vector{String}, params::Parameters; velocities=true)
    # Try to retrieve glacier gdirs if they are available
    filter_missing_glaciers!(rgi_ids, params)
    try
        gdirs::Vector{PyObject} = workflow.init_glacier_directories(rgi_ids)
        filter_missing_glaciers!(gdirs, params)

        # Set different surface topography source if specified
        if params.OGGM.DEM_source != "Default"
            for gdir in gdirs
                tasks.define_glacier_region(gdir, source = params.OGGM.DEM_source)
            end
        end
        
        return gdirs
    catch 
        @warn "Cannot retrieve gdirs from disk!"
        println("Generating gdirs from scratch...")
        # Generate all gdirs if needed
        gdirs::Vector{PyObject} = init_gdirs_scratch(rgi_ids, params; velocities = velocities)
        # Check which gdirs errored in the tasks (useful to filter those gdirs)
        filter_missing_glaciers!(gdirs, params)
        return gdirs
    end
end

"""
    init_gdirs_scratch(rgi_ids)

Initializes Glacier Directories from scratch using OGGM.
"""
function init_gdirs_scratch(rgi_ids::Vector{String}, params::Parameters; velocities=true)::Vector{PyObject}
    # Check if some of the gdirs is missing files
    gdirs::Vector{PyObject} = workflow.init_glacier_directories(rgi_ids, prepro_base_url=params.OGGM.base_url, 
                                                from_prepro_level=2, prepro_border=10,
                                                reset=true, force=true)

    # Set different surface topography source if specified
    if params.OGGM.DEM_source != "Default"
        for gdir in gdirs
            tasks.define_glacier_region(gdir, source = params.OGGM.DEM_source)
        end
    end
    

    if velocities
        list_talks = [
            # tasks.compute_centerlines,
            # tasks.initialize_flowlines,
            # tasks.compute_downstream_line,
            # tasks.catchment_area,
            # tasks.process_dem,
            tasks.gridded_attributes,
            tasks.glacier_masks,
            # tasks.gridded_mb_attributes,
            # tasks.prepare_for_inversion,  # This is a preprocessing task
            # tasks.mass_conservation_inversion,  # This gdirsdoes the actual job
            # tasks.filter_inversion_output,  # This smoothes the thicknesses at the tongue a little
            # tasks.distribute_thickness_per_altitude,
            bedtopo.add_consensus_thickness,   # Use consensus ice thicknesses from Farinotti et al. (2019)
        # tasks.get_topo_predictors,
            millan22.thickness_to_gdir,
            millan22.velocity_to_gdir
        ]
    else
        list_talks = [
            tasks.gridded_attributes,
            tasks.glacier_masks,
            bedtopo.add_consensus_thickness   # Use consensus ice thicknesses from Farinotti et al. (2019)
        ]
    end

    for task in list_talks
        # The order matters!
        workflow.execute_entity_task(task, gdirs)
    end
    GC.gc()

    return gdirs
end

# [Begin] Glathida Utilities
function get_glathida!(gtd_file, gdirs, params; force=false)
    glathida = pd.HDFStore(gtd_file)
    gtd_grids = map(gdir -> get_glathida_glacier(gdir, glathida, force), gdirs)

     # Update missing_glaciers list before removing them
    missing_glaciers = load(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
    for (gtd_grid, gdir) in zip(gtd_grids, gdirs)
        if (length(gtd_grid[gtd_grid .!= 0.0]) == 0) && all(gdir.rgi_id .!= missing_glaciers)
            push!(missing_glaciers, gdir.rgi_id)
            @info "Glacier with all data at 0: $(gdir.rgi_id). Updating list of missing glaciers..."
        end
    end
    jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"); missing_glaciers)

   # Apply deletion to both gtd_grids and gdirs using the same set of indices
    indices_to_remove = findall(x -> length(x[x .!= 0.0]) == 0, gtd_grids)
    deleteat!(gtd_grids, indices_to_remove)
    deleteat!(gdirs, indices_to_remove)

                    
    return gtd_grids, gdirs
end

function get_glathida_glacier(gdir, glathida, force)
    gtd_path = joinpath(gdir.dir, "glathida.h5")
    if isfile(gtd_path) && !force
        gtd_grid = h5read(gtd_path, "gtd_grid")
    else
        df_gtd = glathida[gdir.rgi_id]
        jj, ii = gdir.grid.transform(df_gtd["POINT_LON"], df_gtd["POINT_LAT"], crs=salem.wgs84, nearest=true)

        gtd_grid = zeros((gdir.grid.ny,gdir.grid.nx))
        for (thick, i, j) in zip(df_gtd["THICKNESS"], ii, jj)
            if gtd_grid[i,j] != 0.0
                gtd_grid[i,j] = (gtd_grid[i,j] + thick)/2.0 # average
            else
                gtd_grid[i,j] = thick
            end
        end
        
        # Save file 
        h5open(joinpath(gdir.dir, "glathida.h5"), "w") do file
            write(file, "gtd_grid", gtd_grid)  
        end
    end
    return gtd_grid
end

function get_glathida_path_and_IDs()
    gtd_file = Downloads.download("https://cluster.klima.uni-bremen.de/~oggm/glathida/glathida-v3.1.0/data/TTT_per_rgi_id.h5")
    glathida = pd.HDFStore(gtd_file)
    rgi_ids = glathida.keys()
    rgi_ids = String[id[2:end] for id in rgi_ids]
    return gtd_file, rgi_ids
end
# [End] Glathida Utilities


function filter_missing_glaciers!(gdirs::Vector{PyObject}, params::Parameters)
    task_log::PyObject = global_tasks.compile_task_log(gdirs, 
                                            task_names=["gridded_attributes", "velocity_to_gdir", "thickness_to_gdir"])
                                                        
    task_log.to_csv(joinpath(params.simulation.working_dir, "task_log.csv"))
    glacier_filter = ((task_log.velocity_to_gdir != "SUCCESS").values .&& (task_log.gridded_attributes != "SUCCESS").values
                        .&& (task_log.thickness_to_gdir != "SUCCESS").values)
    glacier_ids = String[]
    for id in task_log.index
        push!(glacier_ids, id)
    end
    missing_glaciers::Vector{String} = glacier_ids[glacier_filter]

    try
        missing_glaciers_old::Vector{String} = load(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
        for missing_glacier in missing_glaciers_old
            if all(missing_glacier .!= missing_glaciers) # if the glacier is already not present, let's add it
                push!(missing_glaciers, missing_glacier)
            end
        end
    catch error
        @warn "$error: No missing_glaciers.jld file available. Skipping..."
    end

    for id in missing_glaciers
        deleteat!(gdirs, findall(x->x.rgi_id==id, gdirs))
    end
    
    # Save missing glaciers in a file
    jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"); missing_glaciers)
    # @warn "Filtering out these glaciers from gdir list: $missing_glaciers"
    
    return missing_glaciers
end

function filter_missing_glaciers!(rgi_ids::Vector{String}, params::Parameters)

    # Check which glaciers we can actually process
    rgi_stats::PyObject = pd.read_csv(utils.file_downloader("https://cluster.klima.uni-bremen.de/~oggm/rgi/rgi62_stats.csv"), index_col=0)
    # rgi_stats = rgi_stats.loc[rgi_ids]

    # if any(rgi_stats.Connect .== 2)
    #     @warn "You have some level 2 glaciers... Removing..."
    #     rgi_ids = [rgi_stats.loc[rgi_stats.Connect .!= 2].index]
    # end

    indices = [rgi_stats.index...]
    for rgi_id in rgi_ids
        if rgi_stats.Connect.values[indices .== rgi_id] == 2
            @warn "Filtering glacier $rgi_id..."
            deleteat!(rgi_ids, rgi_ids .== rgi_id)
        end

    end

    try
        missing_glaciers::Vector{String} = load(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
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

