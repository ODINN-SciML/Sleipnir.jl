export get_rgi_paths

function __init__()

    # Download preprocessed OGGM data
    odinn_path = dirname(prepro_dir)
    if !isdir(odinn_path)
        mkpath(odinn_path)
    end
    if !isdir(prepro_dir)
        @info "Downloading preprocessed data"
        tarGzFile = Downloads.download("https://docs.google.com/uc?export=download&id=1d070a_YqN5aPAONpnzL9hfInv1DA8z3p")
        tar_gz = open(tarGzFile)
        tar = GzipDecompressorStream(tar_gz)
        tempDir = Tar.extract(tar)
        close(tar)
        mv(joinpath(tempDir, "ODINN_prepro"), prepro_dir)
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

function get_rgi_paths()
    rgi_paths = JSON.parsefile(joinpath(prepro_dir, "rgi_paths.json"))
    # Convert Dict{String, Any} to Dict{String, String} and explicitely define type
    # to ensure type stability in the other packages
    rgi_paths::Dict{String, String} = Dict(k => string(v) for (k,v) in pairs(rgi_paths))
    return rgi_paths
end


include("helper_utilities.jl")
