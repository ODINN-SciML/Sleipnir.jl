export netCDF, cfg, utils, workflow, tasks, global_tasks, graphics, bedtopo, millan22, MBsandbox, salem, pd, xr, rioxarray

using Libdl: dlopen

function __init__()

    # Create structural folders if needed
    OGGM_path = joinpath(homedir(), "Python/OGGM_data")
    if !isdir(OGGM_path)
        mkpath(OGGM_path)
    end

    # Avoid issue with dylib files
    try
        if Sys.isapple()
            dlopen(joinpath(root_dir, ".CondaPkg/env/lib/libxml2.dylib"))
            dlopen(joinpath(root_dir, ".CondaPkg/env/lib/libspatialite.8.dylib"))
        elseif Sys.islinux()
            dlopen(joinpath(root_dir, ".CondaPkg/env/lib/libxml2.so.2"))
            load_spatialite()
        else
            error("Unsupported operating system")
        end
    catch e
        @error "Failed to load required libraries" exception=(e, catch_backtrace())
        rethrow(e)
    end
    
    # Load Python packages
    # Only load Python packages if not previously loaded by Sleipnir
    #println("Initializing Python libraries...")
    isassigned(rioxarray) ? nothing : rioxarray[] = pyimport("rioxarray")
    isassigned(netCDF4) ? nothing : netCDF4[] = pyimport("netCDF4")
    isassigned(cfg) ? nothing : cfg[] = pyimport("oggm.cfg")
    isassigned(utils) ? nothing : utils[] = pyimport("oggm.utils")
    isassigned(workflow) ? nothing : workflow[] = pyimport("oggm.workflow")
    isassigned(tasks) ? nothing : tasks[] = pyimport("oggm.tasks")
    isassigned(global_tasks) ? nothing : global_tasks[] = pyimport("oggm.global_tasks")
    isassigned(graphics) ? nothing : graphics[] = pyimport("oggm.graphics")
    isassigned(bedtopo) ? nothing : bedtopo[] = pyimport("oggm.shop.bedtopo")
    isassigned(millan22) ? nothing : millan22[] = pyimport("oggm.shop.millan22")
    isassigned(MBsandbox) ? nothing : MBsandbox[] = pyimport("MBsandbox.mbmod_daily_oneflowline")
    isassigned(salem) ? nothing : salem[] = pyimport("salem")
    isassigned(pd) ? nothing : pd[] = pyimport("pandas")
    isassigned(xr) ? nothing : xr[] = pyimport("xarray")
end

function load_spatialite()
    lib_paths = [
        joinpath(root_dir, ".CondaPkg/env/lib/libspatialite.so"),
        joinpath(root_dir, ".CondaPkg/env/lib/libspatialite.so.2"),
        joinpath(root_dir, ".CondaPkg/env/lib/libspatialite.so.2.12.7")
    ]
    for lib in lib_paths
        if isfile(lib)
            try
                dlopen(lib)
                println("Opened $lib")
            catch e
                println("Failed to load $lib: $e")
            end
        end
    end
end

function clean()
    atexit() do
        run(`$(Base.julia_cmd())`)
    end
    exit()
 end

 function enable_multiprocessing(procs::Int)
    if procs > 0 
        if nprocs() < procs
            @eval begin
            addprocs($procs - nprocs(); exeflags="--project")
            println("Number of cores: ", nprocs())
            println("Number of workers: ", nworkers())
            @everywhere using Sleipnir
            end # @eval
        elseif nprocs() != procs && procs == 1
            @eval begin
            rmprocs(workers(), waitfor=0)
            println("Number of cores: ", nprocs())
            println("Number of workers: ", nworkers())
            end # @eval
        end
    end
    return nworkers()
end



include("helper_utilities.jl")