###############################################################
########  EXTERNAL GEODATA LOADING  ##########################
###############################################################

# ── Hugonnet 2021 geodetic mass balance ──────────────────────

const _HUGONNET_PERIOD = "2000-01-01_2020-01-01"
const _hugonnet_dhdt_cache =
    Ref{Union{Nothing, Dict{String, DhdtData{Sleipnir.Float}}}}(nothing)
const _hugonnet_mb_uncertainty_cache =
    Ref{Union{Nothing, Dict{String, Sleipnir.Float}}}(nothing)

function _default_hugonnet_dhdt_path()
    candidates = (
        joinpath(
            homedir(),
            "OGGM",
            "download_cache",
            "cluster.klima.uni-bremen.de",
            "~oggm",
            "geodetic_ref_mb",
            "hugonnet_2021_ds_rgi60_pergla_rates_10_20_worldwide.csv"),
    )

    for path in candidates
        if isfile(path)
            return path
        end
    end

    return nothing
end

function _parse_hugonnet_period(period::AbstractString)
    bounds = split(period, '_')
    length(bounds) == 2 || return nothing

    start_year = tryparse(Sleipnir.Float, first(split(bounds[1], '-')))
    end_year = tryparse(Sleipnir.Float, first(split(bounds[2], '-')))

    if isnothing(start_year) || isnothing(end_year)
        return nothing
    end

    return (start_year, end_year)
end

function _load_hugonnet_dhdt_cache()
    if !isnothing(_hugonnet_dhdt_cache[])
        return _hugonnet_dhdt_cache[]
    end

    path = _default_hugonnet_dhdt_path()
    if isnothing(path)
        cache = Dict{String, DhdtData{Sleipnir.Float}}()
        _hugonnet_dhdt_cache[] = cache
        return cache
    end

    cache = Dict{String, DhdtData{Sleipnir.Float}}()
    for row in CSV.File(path)
        period_str = string(row.period)
        period_str == _HUGONNET_PERIOD || continue
        ismissing(row.dmdtda) && continue

        period = _parse_hugonnet_period(period_str)
        isnothing(period) && continue

        mb = Sleipnir.Float(row.dmdtda)
        isfinite(mb) || continue

        cache[string(row.rgiid)] = DhdtData(period, mb)
    end

    _hugonnet_dhdt_cache[] = cache
    return cache
end

function _load_hugonnet_mb_uncertainty_cache()
    if !isnothing(_hugonnet_mb_uncertainty_cache[])
        return _hugonnet_mb_uncertainty_cache[]
    end

    path = _default_hugonnet_dhdt_path()
    if isnothing(path)
        cache = Dict{String, Sleipnir.Float}()
        _hugonnet_mb_uncertainty_cache[] = cache
        return cache
    end

    cache = Dict{String, Sleipnir.Float}()
    for row in CSV.File(path)
        period_str = string(row.period)
        period_str == _HUGONNET_PERIOD || continue
        ismissing(row.err_dmdtda) && continue

        period = _parse_hugonnet_period(period_str)
        isnothing(period) && continue

        err_mb = Sleipnir.Float(row.err_dmdtda)
        isfinite(err_mb) || continue

        cache[string(row.rgiid)] = err_mb
    end

    _hugonnet_mb_uncertainty_cache[] = cache
    return cache
end

function _default_hugonnet_dhdt(rgi_id::String)
    return get(_load_hugonnet_dhdt_cache(), rgi_id, nothing)
end

function _default_hugonnet_mb_uncertainty(rgi_id::String)
    return get(_load_hugonnet_mb_uncertainty_cache(), rgi_id, Sleipnir.Float(NaN))
end

# ── GlaThiDa ice thickness ───────────────────────────────────

"""
    get_glathida!(glaciers::Vector{G}, params::Parameters; force=false) where {G <: Glacier2D}

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
function get_glathida!(
        glaciers::Vector{G}, params::Parameters; force = false) where {G <: Glacier2D}
    gtd_grids = pmap(glacier -> get_glathida_glacier(glacier, params, force), glaciers)

    if params.simulation.catch_errors
        # Update missing_glaciers list before removing them
        missing_glaciers = load(joinpath(
            params.simulation.working_dir, "data/missing_glaciers.jld2"))["missing_glaciers"]
        for (gtd_grid, glacier) in zip(gtd_grids, glaciers)
            if (length(gtd_grid[gtd_grid .!= 0.0]) == 0) &&
               all(glacier.rgi_id .!= missing_glaciers)
                push!(missing_glaciers, glacier.rgi_id)
                @info "Glacier with all data at 0: $(glacier.rgi_id). Updating list of missing glaciers..."
            end
        end
        jldsave(joinpath(params.simulation.working_dir, "data/missing_glaciers.jld2");
            missing_glaciers)
    end

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
        for (thick, i, j) in
            zip(glathida["thickness"], glathida["i_grid"], glathida["j_grid"])
            count[i, j] += 1
            gtd_grid[i, j] += thick
        end

        gtd_grid .= ifelse.(count .> 0, gtd_grid ./ count, 0.0)

        # Save file
        h5open(joinpath(rgi_path, "glathida.h5"), "w") do file
            write(file, "gtd_grid", gtd_grid)
        end
    end
    return gtd_grid
end
