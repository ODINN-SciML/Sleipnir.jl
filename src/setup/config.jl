export get_rgi_paths

function __init__()

    # Download preprocessed OGGM data
    odinn_path = dirname(prepro_dir)
    if !isdir(odinn_path)
        mkpath(odinn_path)
    end
    existsAndRedownload = false
    if isdir(prepro_dir)
        daysSinceLastDownload = (Dates.now() - Dates.unix2datetime(mtime(prepro_dir)))/Day(1)
        if daysSinceLastDownload > 7
            # Re-download data if older than one week
            # This is useful especially when the data on the server have been
            # updated and the code needs the new version in order to run
            existsAndRedownload = true
        end
    end
    if (!isdir(prepro_dir)) | existsAndRedownload
        @info "Downloading preprocessed data"
        tarGzFile = Downloads.download("https://drive.usercontent.google.com/download?id=1d070a_YqN5aPAONpnzL9hfInv1DA8z3p&export=download&authuser=1&confirm=t&at=AN8xHooaYmz09qhYCH--PX23RwmA:1753298453038&uuid=f6991b49-5132-4935-8281-5935306e0de1")
        tar_gz = open(tarGzFile)
        tar = GzipDecompressorStream(tar_gz)
        tempDir = Tar.extract(tar)
        close(tar)
        if existsAndRedownload
            rm(prepro_dir, recursive=true)
        end
        mv(joinpath(tempDir, "ODINN_prepro"), prepro_dir)

    end
    csvPath = joinpath(odinn_path, "rgi62_stats.csv")
    if !isfile(csvPath)
        @info "Downloading RGI stats from Bremen cluster"
        csvPathTmp = Downloads.download("https://cluster.klima.uni-bremen.de/~oggm/rgi/rgi62_stats.csv")
        mv(csvPathTmp, csvPath)
    end

end

function clean()
    atexit() do
        run(`$(Base.julia_cmd())`)
    end
    exit()
 end

 function enable_multiprocessing(params::Sleipnir.Parameters)
    procs = params.simulation.workers
    if procs > 0 && params.simulation.multiprocessing
        if nprocs() < procs
            @eval begin
            addprocs($procs - nprocs(); exeflags="--project")
            @info "Number of cores: $(nprocs())"
            @info "Number of workers: $(nworkers())"
            @everywhere using Sleipnir
            end # @eval
        elseif nprocs() != procs && procs == 1
            @eval begin
            rmprocs(workers(), waitfor=0)
            @info "Number of cores: $(nprocs())"
            @info "Number of workers: $(nworkers())"
            end # @eval
        end
    else
        if nprocs()>1 # If the session used to work with multiprocessing but now we want to switch to single processing
            @info "Switching back to single processing"
            @eval begin
            rmprocs(workers(), waitfor=0)
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

function get_rgi_names()
    rgi_names = JSON.parsefile(joinpath(prepro_dir, "rgi_names.json"))
    # Convert Dict{String, Any} to Dict{String, String} and explicitely define type
    # to ensure type stability in the other packages
    rgi_names::Dict{String, String} = Dict(k => string(v) for (k,v) in pairs(rgi_names))
    return rgi_names
end

include("helper_utilities.jl")
