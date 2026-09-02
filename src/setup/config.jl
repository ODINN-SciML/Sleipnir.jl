export get_rgi_paths

function enable_multiprocessing(procs::Int)
    if procs > 0
        if nprocs() < procs
            @eval begin
                project_dir = dirname(Base.active_project())
                # Keep worker logs clean from upstream dependency deprecations.
                addprocs($procs - nprocs(); exeflags = `--project=$(project_dir) --depwarn=no`)
                @info "Number of cores: $(nprocs())"
                @info "Number of workers: $(nworkers())"
                @everywhere using Sleipnir
            end # @eval
        elseif nprocs() != procs && procs == 1
            @eval begin
                rmprocs(workers(), waitfor = 0)
                @info "Number of cores: $(nprocs())"
                @info "Number of workers: $(nworkers())"
            end # @eval
        end
    end
    return nworkers()
end

function get_rgi_paths()
    rgi_paths = JSON.parsefile(joinpath(prepro_dir, "rgi_paths.json"))
    # Convert Dict{String, Any} to Dict{String, String} and explicitely define type
    # to ensure type stability in the other packages
    rgi_paths::Dict{String, String} = Dict(k => string(v) for (k, v) in pairs(rgi_paths))
    return rgi_paths
end

function get_rgi_names()
    rgi_names = JSON.parsefile(joinpath(prepro_dir, "rgi_names.json"))
    # Convert Dict{String, Any} to Dict{String, String} and explicitely define type
    # to ensure type stability in the other packages
    rgi_names::Dict{String, String} = Dict(k => string(v) for (k, v) in pairs(rgi_names))
    return rgi_names
end

include("helper_utilities.jl")
