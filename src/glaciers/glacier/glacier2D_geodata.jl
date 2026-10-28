export glathida_thickness_series

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
The file holds the grid on the native glathida grid, since the `i_grid` and `j_grid` indices of the dataset refer to it. The grid is then averaged by blocks to the size of the glacier when `gridScalingFactor > 1`.
"""
function get_glathida_glacier(glacier::Glacier2D, params::Parameters, force)
    rgi_path = joinpath(prepro_dir(), params.simulation.rgi_paths[glacier.rgi_id])
    gtd_path = joinpath(rgi_path, "glathida.h5")
    if isfile(gtd_path) && !force
        gtd_grid = h5read(gtd_path, "gtd_grid")
    else
        glathida = CSV.File(joinpath(rgi_path, "glathida_data.csv"))
        nx, ny = JSON.parsefile(joinpath(rgi_path, "glacier_grid.json"))["nxny"]
        gtd_grid = zeros(nx, ny)
        count = zeros(nx, ny)
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

    gtd_grid = _coarsen_glathida(gtd_grid, params.simulation.gridScalingFactor)
    @assert size(gtd_grid)==size(glacier.H₀) "The glathida grid of $(glacier.rgi_id) has size $(size(gtd_grid)) which does not match the size $(size(glacier.H₀)) of the ice thickness."
    return gtd_grid
end

# Average a glathida grid by blocks of `n` cells, ignoring the cells without measurement (zero)
function _coarsen_glathida(grid, n)
    n > 1 || return grid
    return block_average_pad_edge_masked(grid, grid .!= 0.0, n; empty_value = 0.0)
end

"""
    _glathida_campaign_index(dates::Vector{Date}, merge_window_days::Integer)

Assign each glathida measurement to a survey campaign.

Two consecutive survey dates belong to the same campaign when they are at most
`merge_window_days` apart, which groups the multi day field campaigns that are reported as
separate dates in the raw dataset. The rule is applied date by date, so a series of surveys
that are each less than `merge_window_days` apart is one campaign.

# Arguments

  - `dates::Vector{Date}`: Date of each measurement.
  - `merge_window_days::Integer`: Maximum number of days between two surveys of a campaign.

# Returns

  - `index`: A vector of campaign indices, one per measurement, increasing with time.
"""
function _glathida_campaign_index(dates::Vector{Date}, merge_window_days::Integer)
    isempty(dates) && throw(ArgumentError("There are no dated glathida measurements."))
    unique_dates = sort(unique(dates))
    campaign_of_date = Dict{Date, Int}()
    campaign = 1
    campaign_of_date[unique_dates[1]] = campaign
    for k in 2:length(unique_dates)
        if Dates.value(unique_dates[k] - unique_dates[k - 1]) > merge_window_days
            campaign += 1
        end
        campaign_of_date[unique_dates[k]] = campaign
    end
    return [campaign_of_date[d] for d in dates]
end

"""
    get_glathida_campaigns(glacier::Glacier2D, params::Parameters, force; merge_window_days::Integer = 30)

Retrieve or generate the per campaign glathida grids for a given glacier.

[`get_glathida_glacier`](@ref) averages all the measurements of a cell into one dateless
grid. This function does not: it returns one grid for each survey campaign, so that each
campaign is a separate observation of a transient inversion. A cell measured in two campaigns
has two values, one in each grid. Only the measurements of the same cell and the same
campaign are averaged. Surveys at most `merge_window_days` apart are one campaign.

# Arguments

  - `glacier::Glacier2D`: The glacier object for which the glathida campaigns are to be retrieved or generated.
  - `params::Parameters`: The parameters object containing simulation settings.
  - `force`: A boolean flag indicating whether to force regeneration of the campaigns even if they already exist.
  - `merge_window_days::Integer`: Maximum number of days between two surveys of a campaign.

# Returns

  - `t`: A vector of decimal years, one per campaign, sorted chronologically.
  - `gtd_grids`: A vector of 2D arrays on the native glathida grid, with zeros where no measurement is available.
  - `counts`: The number of measurements in each cell, one 2D array per campaign.
  - `sumsqs`: The sum of the squared thicknesses in each cell, one 2D array per campaign. With
    `counts`, it gives the spread of the measurements inside a cell.

# Description

This function checks if the glathida campaign file (`glathida_campaigns_<merge_window_days>d.h5`)
exists in the specified path. If the file exists and `force` is `false`, it reads the
campaigns from the file. Otherwise, it reads the glacier thickness data from a CSV file
(`glathida_data.csv`), groups the measurements by campaign, computes the average thickness
for each grid cell of each campaign, and saves the result to this HDF5 file.
The grids are built on the native glathida grid since the `i_grid` and `j_grid` indices of
the dataset refer to it, irrespective of `gridScalingFactor`.
"""
function get_glathida_campaigns(
        glacier::Glacier2D,
        params::Parameters,
        force;
        merge_window_days::Integer = 30
)
    rgi_path = joinpath(prepro_dir(), params.simulation.rgi_paths[glacier.rgi_id])
    gtd_path = joinpath(rgi_path, "glathida_campaigns_$(merge_window_days)d.h5")
    # A file written before the counts were saved is built again
    cached = isfile(gtd_path) && !force && h5open(f -> haskey(f, "counts"), gtd_path)
    if cached
        t = h5read(gtd_path, "t")
        gtd_grids = h5read(gtd_path, "gtd_grids")
        counts = h5read(gtd_path, "counts")
        sumsqs = h5read(gtd_path, "sumsqs")
    else
        glathida = CSV.File(joinpath(rgi_path, "glathida_data.csv"))
        nx, ny = JSON.parsefile(joinpath(rgi_path, "glacier_grid.json"))["nxny"]

        dates = Date.(glathida["date"])
        campaign = _glathida_campaign_index(dates, merge_window_days)
        ncampaigns = maximum(campaign)

        gtd_grids = zeros(nx, ny, ncampaigns)
        counts = zeros(nx, ny, ncampaigns)
        sumsqs = zeros(nx, ny, ncampaigns)
        for (c, thick, i, j) in
            zip(campaign, glathida["thickness"], glathida["i_grid"], glathida["j_grid"])
            counts[i, j, c] += 1
            gtd_grids[i, j, c] += thick
            sumsqs[i, j, c] += thick^2
        end

        gtd_grids .= ifelse.(counts .> 0, gtd_grids ./ counts, 0.0)

        # Date of a campaign: the mean date of its measurements
        t = [mean(datetime_to_floatyear.(DateTime.(dates[campaign .== c])))
             for c in 1:ncampaigns]

        # Save file
        h5open(gtd_path, "w") do file
            write(file, "t", t)
            write(file, "gtd_grids", gtd_grids)
            write(file, "counts", counts)
            write(file, "sumsqs", sumsqs)
        end
    end
    split(a) = [a[:, :, c] for c in axes(a, 3)]
    return t, split(gtd_grids), split(counts), split(sumsqs)
end

"""
    glathida_cell_σ(H, grid, count, sumsq, n; σ_picking = 5.0, σ_velocity = 0.03,
        representation_factor = 1.5, σ_spread = 12.0, min_points = 3)

Uncertainty of each cell of a glathida campaign averaged by blocks of `n` cells, in m:

``σ^2 = σ_{picking}^2 + (σ_{velocity} H)^2 + (f s)^2``

The first two parts are the error of a radar measurement. They come from the picking of the
bed reflection (`σ_picking`, in m) and from the radio-wave velocity (`σ_velocity`, as a
fraction of the thickness). The two errors are independent, so they add in quadrature
(Lapazaran et al., 2016, Eq. 5). The defaults are those of Grab et al. (2021): ±5 m for the
picking on average (±3 m for a well defined reflection, ±7 m for a poorly defined one) and
0.169 ± 0.005 m ns⁻¹ for the velocity in temperate ice, which is 3%. Lapazaran et al. (2016)
give 1 to 5% for the velocity, and a timing error that depends on the frequency of the radar
(4.2 m at 20 MHz, 0.42 m at 200 MHz).

The last part is the representation error: the model cell is the mean over the whole cell,
while a survey line only crosses it. `s` is the standard deviation of the measurements inside
the cell, and the factor `f = representation_factor` is there because one line sees only part
of the cell. Its default was measured on Aletsch with two independent campaigns, at 137, 274
and 548 m. Cells with fewer than `min_points` measurements take the median `s` of the
campaign, or `σ_spread` (in m) if no cell of the campaign has `min_points` measurements.

`H` is the averaged grid. `grid` (mean of each cell), `count` and `sumsq` are on the native
grid. The cells without measurement get 0.

# References

  - Lapazaran, J. J., Otero, J., Martín-Español, A. and Navarro, F. J. (2016). On the errors
    involved in ice-thickness estimates I: ground-penetrating radar measurement errors.
    Journal of Glaciology, 62(236), 1008-1020.
  - Grab, M. et al. (2021). Ice thickness distribution of all Swiss glaciers based on extended
    ground-penetrating radar data and glaciological modeling. Journal of Glaciology, 67(266),
    1074-1092.
"""
function glathida_cell_σ(H, grid, count, sumsq, n; σ_picking::Real = 5.0,
        σ_velocity::Real = 0.03, representation_factor::Real = 1.5, σ_spread::Real = 12.0,
        min_points::Integer = 3)
    min_points >= 2 ||
        throw(ArgumentError("min_points must be at least 2 to have a spread."))
    coarse(a) = n > 1 ? block_average_pad_edge(a, n) : a
    # The block average divides the sums and the count by the same number, so the mean and
    # the mean square of the measurements of a cell are their ratios
    N = coarse(count)
    mean_H = coarse(grid .* count) ./ N
    mean_H² = coarse(sumsq) ./ N
    npoints = round.(Int, N .* n^2)
    observed = (H .!= 0.0) .& (N .> 0)
    s = fill(NaN, size(H))
    many = observed .& (npoints .>= min_points)
    s[many] .= sqrt.(max.(mean_H²[many] .- mean_H[many] .^ 2, 0.0) .*
                     npoints[many] ./ (npoints[many] .- 1))
    # With no cell of `min_points` measurements the spread is not known, and the radar error
    # alone (5 to 10 m) would be too small: 12 m gives about 19 m for a cell of 140 m of ice,
    # close to the 18 m measured on Aletsch
    s[observed .& .!many] .= any(many) ? median(s[many]) : σ_spread
    σ = zeros(size(H))
    σ[observed] .= sqrt.(σ_picking^2 .+ (σ_velocity .* H[observed]) .^ 2 .+
                         (representation_factor .*
                          ifelse.(isnan.(s[observed]), 0.0, s[observed])) .^ 2)
    return σ
end

"""
    glathida_thickness_series(glacier::Glacier2D, params::Parameters; force = false, merge_window_days::Integer = 30, kwargs...)

Build the ice thickness time series of the glathida survey campaigns of a glacier.

The returned [`ThicknessData`](@ref) can be attached to a glacier so that each campaign
constrains the ice thickness at its own date during a transient inversion. Unobserved cells
are left at zero, which is the convention `is_in_glacier` relies on to select the observed
pixels. Since the resulting grids are sparse, the accompanying loss must use a zero distance
to the border, otherwise the mask is eroded to nothing.

# Arguments

  - `glacier::Glacier2D`: The glacier object for which the thickness series is built.
  - `params::Parameters`: The parameters object containing simulation settings.
  - `force`: A boolean flag indicating whether to force regeneration of the campaigns even if they already exist.
  - `merge_window_days::Integer`: Maximum number of days between two surveys of a campaign.
  - `kwargs...`: The keywords of [`glathida_cell_σ`](@ref) (`σ_picking`, `σ_velocity`,
    `representation_factor`, `σ_spread` and `min_points`), to change the uncertainty of the cells.

# Returns

  - `thicknessData`: A `ThicknessData` whose entries are the glathida campaigns, sorted
    chronologically, with the uncertainty of each cell from [`glathida_cell_σ`](@ref).
"""
function glathida_thickness_series(
        glacier::Glacier2D,
        params::Parameters;
        force = false,
        merge_window_days::Integer = 30,
        kwargs...
)
    t, gtd_grids,
    counts,
    sumsqs = get_glathida_campaigns(
        glacier, params, force; merge_window_days = merge_window_days)

    n = params.simulation.gridScalingFactor
    H = _coarsen_glathida.(gtd_grids, n)
    @assert all(size.(H) .== Ref(size(glacier.H₀))) "The glathida grids of $(glacier.rgi_id) have size $(size(first(H))) which does not match the size $(size(glacier.H₀)) of the ice thickness."
    σ = [glathida_cell_σ(H[c], gtd_grids[c], counts[c], sumsqs[c], n; kwargs...)
         for c in eachindex(H)]

    return ThicknessData(t = t, H = H, σ = σ)
end
