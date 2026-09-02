###############################################################
########  EXTERNAL GEODATA LOADING  ##########################
###############################################################

# ── Hugonnet 2021 geodetic mass balance ──────────────────────

const _HUGONNET_PERIOD = (Date(2000, 1, 1), Date(2020, 1, 1))
# Decimal-year bounds of the Hugonnet calibration window. The raw climate series must
# span this period when the mass balance model is calibrated against geodetic data.
const HUGONNET_CLIMATE_PERIOD = (2000.0, 2020.0)
const _hugonnet_df_cache = Ref{Union{
    Nothing, DataFrame}}(nothing)

function _load_hugonnet_dataset()
    path = joinpath(artifact"hugonnet21_dataset", "hugonnet21_dataset",
        "hugonnet_2021_ds_rgi60_pergla_rates_10_20_worldwide_filled.csv")
    df = CSV.read(path, DataFrame)
    df.period_start = Date.(first.(split.(df.period, "_")))
    df.period_end = Date.(last.(split.(df.period, "_")))
    df = subset(df, :period_start => ByRow(==(_HUGONNET_PERIOD[1])),
        :period_end => ByRow(==(_HUGONNET_PERIOD[2])))
    df = df[.!ismissing.(df.dmdtda), :]
    df = df[.!ismissing.(df.err_dmdtda), :]
    return df
end

function _build_dhdt(df::DataFrame, rgi_id::String)
    if rgi_id in df.rgiid
        sel_df = df[df.rgiid .== rgi_id, :]
        @assert size(sel_df, 1) == 1 "Found multiple entries for $(rgiid) in the Hugonnet21 dataset for the geodetic period $(_HUGONNET_PERIOD)."
        period = (Sleipnir.Float(year(sel_df[1, :period_start])),
            Sleipnir.Float(year(sel_df[1, :period_end])))
        mb = Sleipnir.Float(sel_df[1, :dmdtda])
        err_mb = Sleipnir.Float(sel_df[1, :err_dmdtda])
        return DhdtData(period, mb, err_mb)
    else
        return nothing
    end
end

function _load_hugonnet_df_cache()
    isnothing(_hugonnet_df_cache[]) || return _hugonnet_df_cache[]
    df = _load_hugonnet_dataset()
    _hugonnet_df_cache[] = df
    return df
end

function _default_hugonnet_dhdt(rgi_id::String)
    return _build_dhdt(_load_hugonnet_df_cache(), rgi_id)
end
function _default_hugonnet_dhdt(rgi_ids::Vector{String})
    map(rgi_ids) do rgi_id
        _build_dhdt(_load_hugonnet_df_cache(), rgi_id)
    end
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
