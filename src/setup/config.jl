
function __init__()

    # Create structural folders if needed
    OGGM_path = joinpath(homedir(), "Python/OGGM_data")
    if !isdir(OGGM_path)
        mkpath(OGGM_path)
    end

    # # Avoid issue with dylib files
    # try
    #     if Sys.isapple()
    #         dlopen(joinpath(root_dir, ".CondaPkg/env/lib/libxml2.dylib"))
    #         dlopen(joinpath(root_dir, ".CondaPkg/env/lib/libspatialite.8.dylib"))
    #     elseif Sys.islinux()
    #         load_libxml()
    #         load_spatialite()
    #     else
    #         error("Unsupported operating system")
    #     end
    # catch e
    #     @error "Failed to load required libraries" exception=(e, catch_backtrace())
    #     rethrow(e)
    # end
end

# function load_libxml()
#     lib_dir = joinpath(root_dir, ".CondaPkg/env/lib")
    
#     # Find all libspatialite files in the directory
#     lib_files = filter(f -> startswith(f, "libxml") && (endswith(f, ".so") || contains(f, ".so.")), readdir(lib_dir))
    
#     if isempty(lib_files)
#         println("No libxml files found in $lib_dir")
#         return
#     end
    
#     for lib_file in lib_files
#         lib_path = joinpath(lib_dir, lib_file)
#         try
#             dlopen(lib_path)
#             println("Opened $lib_path")
#         catch e
#             println("Failed to load $lib_path: $e")
#         end
#     end
# end

# function load_spatialite()
#     lib_dir = joinpath(root_dir, ".CondaPkg/env/lib")
    
#     # Find all libspatialite files in the directory
#     lib_files = filter(f -> startswith(f, "libspatialite") && (endswith(f, ".so") || contains(f, ".so.")), readdir(lib_dir))
    
#     if isempty(lib_files)
#         println("No libspatialite files found in $lib_dir")
#         return
#     end
    
#     for lib_file in lib_files
#         lib_path = joinpath(lib_dir, lib_file)
#         try
#             dlopen(lib_path)
#             println("Opened $lib_path")
#         catch e
#             println("Failed to load $lib_path: $e")
#         end
#     end
# end

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


# function load_lib(libname::String)
    
#     # Find all libspatialite files in the directory

#     # Find way to pass this path
#     lib_dir = joinpath(dirname(dirname(read(`which python`, String)[1:end-1])), "lib")
#     # lib_dir = "/usr/local/Caskroom/miniforge/base/envs/oggm_env_20240917_ssl/lib"

#     if Sys.isapple()
#         lib_files = filter(f -> startswith(f, libname) && (endswith(f, ".dylib") || contains(f, ".dylib.")), readdir(lib_dir))
#     elseif Sys.islinux()
#         lib_files = filter(f -> startswith(f, libname) && (endswith(f, ".so") || contains(f, ".so.")), readdir(lib_dir))
#     else
#         error("Unsupported operating system")
#     end

#     if isempty(lib_files)
#         println("No libxml files found in $lib_dir")
#         return
#     end
    
#     for lib_file in lib_files
#         lib_path = joinpath(lib_dir, lib_file)
#         try
#             dlopen(lib_path)
#             println("Opened $lib_path")
#         catch e
#             println("Failed to load $lib_path: $e")
#         end
#     end
# end

include("helper_utilities.jl")
