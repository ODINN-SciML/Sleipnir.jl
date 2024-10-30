export netCDF4, cfg, utils, workflow, tasks, global_tasks, graphics, bedtopo, millan22, MBsandbox, salem, pd, xr, rioxarray
using Libdl: dlopen

function __init__()

    # Create structural folders if needed
    OGGM_path = joinpath(homedir(), "Python/OGGM_data")
    if !isdir(OGGM_path)
        mkpath(OGGM_path)
    end

    # Avoid issue with dylib files
    # try
        load_lib("libxml")
        load_lib("libspatialite")
        load_lib("libcrypto")
    # catch e
        @error "Failed to load required libraries" exception=(e, catch_backtrace())
        rethrow(e)
    # end

    # Load Python packages
    # try
    #     # Only load Python packages if not previously loaded by Sleipnir
        if cfg == PyNULL() && workflow == PyNULL() && utils == PyNULL() && MBsandbox == PyNULL() 
            println("Initializing Python libraries...")
            copy!(netCDF4, pyimport("netCDF4"))
            copy!(cfg, pyimport("oggm.cfg"))
            copy!(utils, pyimport("oggm.utils"))
            copy!(workflow, pyimport("oggm.workflow"))
            copy!(tasks, pyimport("oggm.tasks"))
            copy!(global_tasks, pyimport("oggm.global_tasks"))
            copy!(graphics, pyimport("oggm.graphics"))
            copy!(bedtopo, pyimport("oggm.shop.bedtopo"))
            copy!(millan22, pyimport("oggm.shop.millan22"))
            copy!(MBsandbox, pyimport("MBsandbox.mbmod_daily_oneflowline"))
            copy!(salem, pyimport("salem"))
            copy!(pd, pyimport("pandas"))
            copy!(xr, pyimport("xarray"))
            copy!(rioxarray, pyimport("rioxarray"))
        end
    # catch e
    #     @warn "It looks like you have not installed and/or activated the virtual Python environment. \n 
    #     Please follow the guidelines in: https://github.com/ODINN-SciML/ODINN.jl#readme"
    #     @warn exception=(e, catch_backtrace())
    # end

    if isempty(lib_files)
        println("No libxml files found in $lib_dir")
        return
    end
    
    for lib_file in lib_files
        lib_path = joinpath(lib_dir, lib_file)
        try
            dlopen(lib_path)
            println("Opened $lib_path")
        catch e
            println("Failed to load $lib_path: $e")
        end
    end
end

function load_spatialite()
    lib_dir = joinpath(root_dir, ".CondaPkg/env/lib")
    # @show lib_dir
    
    # Find all libspatialite files in the directory
    if Sys.isapple()
        lib_files = filter(f -> startswith(f, "libspatialite") && (endswith(f, ".dylib") || contains(f, ".dylib.")), readdir(lib_dir))
    elseif Sys.islinux()
        lib_files = filter(f -> startswith(f, "libspatialite") && (endswith(f, ".so") || contains(f, ".so.")), readdir(lib_dir))
    else
        error("Unsupported operating system")
    end

    if isempty(lib_files)
        println("No libspatialite files found in $lib_dir")
        return
    end
    
    for lib_file in lib_files
        lib_path = joinpath(lib_dir, lib_file)
        try
            dlopen(lib_path)
            println("Opened $lib_path")
        catch e
            println("Failed to load $lib_path: $e")
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

function filter_existing_paths(paths::Vector{String})
    # Use `filter` to retain only the paths that exist
    existing_paths = filter(ispath, paths)
    return existing_paths
end


function load_lib(libname::String)
    
    # Find all libspatialite files in the directory

    # Find way to pass this path
    lib_dir = joinpath(dirname(dirname(read(`which python`, String)[1:end-1])), "lib")
    # lib_dir = "/usr/local/Caskroom/miniforge/base/envs/oggm_env_20240917_ssl/lib"

    if Sys.isapple()
        lib_files = filter(f -> startswith(f, libname) && (endswith(f, ".dylib") || contains(f, ".dylib.")), readdir(lib_dir))
    elseif Sys.islinux()
        lib_files = filter(f -> startswith(f, libname) && (endswith(f, ".so") || contains(f, ".so.")), readdir(lib_dir))
    else
        error("Unsupported operating system")
    end

    if isempty(lib_files)
        println("No libxml files found in $lib_dir")
        return
    end
    
    for lib_file in lib_files
        lib_path = joinpath(lib_dir, lib_file)
        try
            dlopen(lib_path)
            println("Opened $lib_path")
        catch e
            println("Failed to load $lib_path: $e")
        end
    end
end

include("helper_utilities.jl")
