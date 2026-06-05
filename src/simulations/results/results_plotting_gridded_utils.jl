###############################################################################
# Gridded-data plotting utilities (CairoMakie)
#
# Layout:
#   1. Numeric helpers:   lower_bound, finite_extrema, symmetric_finite_colorrange
#   2. Internal helpers:  _resolve_gridded_snapshot, _dem_defaults
#   3. Core renderer:     _plot_gridded_data_core
#   4. Public API:        plot_gridded_data, accumulate_gridded_data,
#                         plot_cumulative_gridded_data, plot_cumulative_mb,
#                         plot_glacier_dem, save_figure
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 1. Numeric helpers
# ─────────────────────────────────────────────────────────────────────────────

function lower_bound(M::Matrix{<: AbstractFloat})
    positive_values = M[(!isnan).(M) .& (M .> 0)]
    if !isempty(positive_values)
        return minimum(positive_values)
    end

    finite_values = filter(isfinite, vec(M))
    @assert !isempty(finite_values) "There are no finite values to plot."
    return minimum(finite_values)
end

function finite_extrema(M::Matrix{<: AbstractFloat})
    values = filter(isfinite, vec(M))
    @assert !isempty(values) "There are no finite values to plot."
    return extrema(values)
end

function symmetric_finite_colorrange(M::Matrix{<: AbstractFloat})
    values = filter(isfinite, vec(M))
    @assert !isempty(values) "There are no finite values to plot."

    max_abs = maximum(abs, values)
    max_abs = max(max_abs, eps(Float64))
    return (-max_abs, max_abs)
end

# ─────────────────────────────────────────────────────────────────────────────
# 2. Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

function _resolve_gridded_snapshot(
        gridded_data::Union{Vector{Matrix{F}}, Matrix{F}},
        timeIdx::Union{Nothing, Int64};
        data_name::String = "Data"
) where {F <: AbstractFloat}
    if gridded_data isa AbstractVector
        @assert !isempty(gridded_data) "$(data_name) is an empty vector"
        @assert isnothing(timeIdx) || length(gridded_data) >= timeIdx "The provided index=$(timeIdx) is greater than the size of the vector which is $(length(gridded_data))"
        return isnothing(timeIdx) ? gridded_data[end] : gridded_data[timeIdx]
    end
    return gridded_data
end

function _dem_defaults(results::Results)
    return (
        title = "Glacier DEM",
        colorbar_label = "Surface elevation (m a.s.l.)"
    )
end

function _dem_defaults(glacier::Glacier2D)
    plot_title = isempty(glacier.name) ? "Glacier DEM" : "Glacier DEM ($(glacier.name))"
    return (
        title = plot_title,
        colorbar_label = "Surface elevation (m a.s.l.)"
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# 3. Core renderer
# ─────────────────────────────────────────────────────────────────────────────

function _plot_gridded_data_core(
        gridded_data::Union{Vector{Matrix{F}}, Matrix{F}};
        lon::F,
        lat::F,
        x::Vector{F},
        y::Vector{F},
        rgi_id::String,
        Δx::F,
        mask::BitMatrix,
        scale_text_size::Union{Nothing, Float64} = nothing,
        timeIdx::Union{Nothing, Int64} = nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        plotContour::Bool = false,
        colormap = :cool,
        colorrange::Union{Nothing, Tuple{<: Real, <: Real}} = nothing,
        colorscale = nothing,
        logPlot = false,
        title::Union{Nothing, String} = nothing,
        colorbar_label::Union{Nothing, String} = nothing
) where {F <: AbstractFloat}
    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)

    nx, ny = size(mask)
    ctr = plotContour ? Contour.contour(collect(1:nx), 1+ny .- collect(1:ny), mask, 0.5) :
          nothing

    current_data = _resolve_gridded_snapshot(gridded_data, timeIdx; data_name = "Data")

    @assert size(current_data) == size(mask) "The data and mask must have the same size."

    current_data_masked = copy(current_data)
    current_data_masked[.!mask] .= NaN

    if isnothing(colorrange)
        if logPlot
            global_min = lower_bound(current_data_masked)
            global_max = maximum(replace(current_data_masked, NaN => 0.0))
        else
            global_min, global_max = finite_extrema(current_data_masked)
        end
    else
        global_min, global_max = colorrange
    end

    figKwargs[:layout] = GridLayout(2, 2)
    fig = Figure(; figKwargs...)

    ax_row = 1
    ax_col = 1
    ax = Axis(fig[ax_row, ax_col], aspect = DataAspect())
    data = copy(current_data)

    nx, ny = size(data)

    data[.!mask] .= NaN

    hm = heatmap!(ax, reverseForHeatmap(data, x, y), colormap = colormap,
        colorrange = (logPlot ? global_min*0.85 : global_min, global_max),
        colorscale = isnothing(colorscale) ? (logPlot ? log10 : identity) : colorscale)
    cb_kwargs = isnothing(colorbar_label) ? NamedTuple() : (; label = colorbar_label)
    cb = Colorbar(fig[ax_row, ax_col + 1], hm; cb_kwargs...)
    Observables.connect!(cb.height, @lift CairoMakie.Fixed($(viewport(ax.scene)).widths[2]))

    if plotContour
        _overlay_contour!(ax, ctr)
    end

    _decorate_geo_axis!(ax, nx, ny, lon, lat, Δx;
        scale_text_size = scale_text_size)

    fig_title = isnothing(title) ? rgi_id : "$title — $rgi_id"
    fig[0, :] = Label(fig, fig_title, fontsize = 14, font = :bold)
    resize_to_layout!(fig)
    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
# 4. Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_gridded_data(
        gridded_data::Union{Vector{Matrix{F}}, Matrix{F}},
        results::Results;
        scale_text_size::Union{Nothing,Float64}=nothing,
        timeIdx::Union{Nothing,Int64}=nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        plotContour::Bool=false,
        colormap = :cool,
        logPlot = false,
    ) where {F <: AbstractFloat}

Plot a gridded matrix (or a time series of matrices) as a heatmap using metadata from results.

# Arguments

  - `gridded_data::Union{Vector{Matrix{F}}, Matrix{F}}`: Single snapshot or time series (defaults to last timestep).
  - `results::Results`: Supplies lon, lat, x, y, rgi_id, Δx and H (mask).
  - `scale_text_size`, `figsize`, `colormap`: Optional plotting params.
  - `timeIdx::Union{Nothing,Int64}`: Select timestep when `gridded_data` is a vector.
  - `plotContour::Bool`: overlay glacier-mask contour from results.H.
  - `logPlot::Bool`: Use log10 colorscale (positive non-NaN values determine range).

# Behavior

  - Masks out cells where `results.H[begin] .<= 0` (set to NaN).
  - Adds colorbar, central lon/lat tick, and a Δx-based scale bar in km.
  - If `plotContour`, draws mask boundary lines.
  - Returns a `CairoMakie.Figure`.

# Errors

  - Asserts gridded_data is non-empty and timeIdx (if provided) is in range.
"""
function plot_gridded_data(
        gridded_data::Union{Vector{Matrix{F}}, Matrix{F}},
        results::Results;
        scale_text_size::Union{Nothing, Float64} = nothing,
        timeIdx::Union{Nothing, Int64} = nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        plotContour::Bool = false,
        colormap = :cool,
        colorrange::Union{Nothing, Tuple{<: Real, <: Real}} = nothing,
        colorscale = nothing,
        logPlot = false,
        title::Union{Nothing, String} = nothing,
        colorbar_label::Union{Nothing, String} = nothing
) where {F <: AbstractFloat}
    mask = results.H[begin] .> 0.0
    return _plot_gridded_data_core(gridded_data;
        lon = results.lon,
        lat = results.lat,
        x = results.x,
        y = results.y,
        rgi_id = results.rgi_id,
        Δx = results.Δx,
        mask = mask,
        scale_text_size = scale_text_size,
        timeIdx = timeIdx,
        figsize = figsize,
        plotContour = plotContour,
        colormap = colormap,
        colorrange = colorrange,
        colorscale = colorscale,
        logPlot = logPlot,
        title = title,
        colorbar_label = colorbar_label)
end

"""
    accumulate_gridded_data(
        gridded_data::Vector{Matrix{F}};
        weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    ) where {F <: AbstractFloat}

Accumulate a time series of gridded matrices into a single matrix.

# Arguments

  - `gridded_data::Vector{Matrix{F}}`: Sequence of gridded fields to accumulate.
  - `weights::Union{Nothing,AbstractVector{<:Real}}`: Optional per-step weights. If
    provided, weighted accumulation is performed as `sum(weights[i] * gridded_data[i])`.

# Returns

  - `Matrix{F}`: The accumulated matrix.
"""
function accumulate_gridded_data(
        gridded_data::Vector{Matrix{F}};
        weights::Union{Nothing, AbstractVector{<:Real}} = nothing
) where {F <: AbstractFloat}
    @assert !isempty(gridded_data) "There is no data to accumulate."

    cumulative = zeros(F, size(gridded_data[begin]))
    if isnothing(weights)
        cumulative = sum(gridded_data)
    else
        @assert length(weights) == length(gridded_data) "Length of weights must match gridded_data."
        for i in eachindex(gridded_data)
            cumulative .+= weights[i] .* gridded_data[i]
        end
    end

    return cumulative
end

"""
    plot_cumulative_gridded_data(
        gridded_data::Vector{Matrix{F}},
        results::Results;
        weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
        kwargs...
    ) where {F <: AbstractFloat}

Plot the cumulative field of a time series of gridded matrices using `plot_gridded_data`.
This is a thin utility wrapper to avoid duplicating plotting code.

# Arguments

  - `gridded_data::Vector{Matrix{F}}`: Sequence of gridded fields to accumulate.
  - `results::Results`: Results object with glacier metadata for plotting.
  - `weights::Union{Nothing,AbstractVector{<:Real}}`: Optional per-step weights.
  - `kwargs...`: Additional keyword arguments forwarded to `plot_gridded_data`.

# Returns

  - `Figure`: Cumulative field figure.
"""
function plot_cumulative_gridded_data(
        gridded_data::Vector{Matrix{F}},
        results::Results;
        weights::Union{Nothing, AbstractVector{<:Real}} = nothing,
        kwargs...
) where {F <: AbstractFloat}
    cumulative = accumulate_gridded_data(gridded_data; weights = weights)
    return plot_gridded_data(cumulative, results; kwargs...)
end

"""
    plot_cumulative_mb(results::Results; kwargs...)

Plot a spatial map of the **cumulative** surface mass balance accumulated over the
simulation period from the per-callback MB fields stored in `results.MB`.

Each entry of `results.MB` is the gridded mass-balance increment produced at one MB
callback (one per MB time step). The function sums these increments cell-by-cell
(via `accumulate_gridded_data`) to obtain the total mass balance over the period
`results.tspan = (t0, t1)`.

# Modes (`annual_MB`)

  - `annual_MB = false` (default): plots the raw cumulative MB over the whole period,
    i.e. `Σᵢ MBᵢ`, in `m w.e.`.
  - `annual_MB = true`: divides every increment by the period length `T = t1 - t0`
    (in years) before summing, i.e. `Σᵢ (MBᵢ / T) = (Σᵢ MBᵢ) / T`. Since `Σᵢ MBᵢ` is the
    total over `T` years, this is the **mean annual mass-balance rate** (the
    "annually-averaged equivalent"), in `m w.e. yr⁻¹`.

Returns `nothing` (with a warning) when `results.MB` is empty — e.g. when the
simulation was run without `use_MB = true`.

# Arguments

  - `results::Results`: results carrying the per-callback MB maps in `results.MB`.

# Keyword arguments

  - `title::String`: figure title prefix (the period `(t0–t1)` is appended automatically).
  - `colorbar_label::Union{Nothing,String}`: colorbar label; defaults to `m w.e.`
    (or `m w.e. yr⁻¹` when `annual_MB = true`).
  - `annual_MB::Bool = false`: switch between cumulative total and mean annual rate.
  - `colormap`: diverging colormap (red→white→blue by default).
  - `kwargs...`: forwarded to [`plot_gridded_data`](@ref).

# Returns

  - `Figure`, or `nothing` if there is no MB history.
"""
function plot_cumulative_mb(results::Results;
        title::String = "Cumulative Mass Balance",
        colorbar_label::Union{Nothing, String} = nothing,
        annual_MB::Bool = false,
        colormap = CairoMakie.cgrad([:indianred3, :white, :steelblue3], [0.0, 0.5, 1.0]),
        kwargs...)
    if isempty(results.MB) || isempty(results.MB[begin])
        @warn "No mass balance callback history in results; skipping cumulative MB plot. " *
              "Make sure the simulation was run with use_MB=true."
        return nothing
    end
    mb_data = results.MB
    t0, t1 = results.tspan
    period_str = "$(round(Int, t0))–$(round(Int, t1))"
    full_title = "$title ($period_str)"
    if annual_MB
        duration = t1 - t0
        if duration > 0
            mb_data = [m ./ duration for m in mb_data]
        end
        label = isnothing(colorbar_label) ? "m w.e. yr⁻¹" : colorbar_label
    else
        label = isnothing(colorbar_label) ? "m w.e." : colorbar_label
    end
    cumulative_mb = accumulate_gridded_data(mb_data)
    mb_mask = results.H[begin] .> 0.0
    cumulative_mb_masked = copy(cumulative_mb)
    cumulative_mb_masked[.!mb_mask] .= NaN
    mb_colorrange = symmetric_finite_colorrange(cumulative_mb_masked)
    return plot_gridded_data(cumulative_mb, results;
        title = full_title,
        colorbar_label = label,
        colormap = colormap,
        colorrange = mb_colorrange,
        kwargs...)
end

"""
    plot_glacier_dem(glacier_or_results; kwargs...)

Plot the glacier DEM (surface elevation field `S`) with a terrain colormap,
geographic coordinate labels, colorbar, and glacier contour overlay.

# Arguments

    - `glacier_or_results`: Either `Results` or `Glacier2D`.

# Keyword Arguments

    - `title::String`: Figure title prefix.
    - `colorbar_label::String`: Label for the colorbar.
    - `plotContour::Bool`: Whether to overlay glacier contour (default `true`).
    - `colormap`: Makie colormap symbol (default `:terrain`).
    - `kwargs...`: Additional keyword arguments forwarded to `plot_gridded_data`.

# Returns

    - `Figure`: DEM figure.
"""
function plot_glacier_dem(results::Results;
        title::Union{Nothing, String} = nothing,
        colorbar_label::Union{Nothing, String} = nothing,
        plotContour::Bool = true,
        colormap = :terrain,
        kwargs...)
    defaults = _dem_defaults(results)

    return plot_gridded_data(results.S, results;
        title = isnothing(title) ? defaults.title : title,
        colorbar_label = isnothing(colorbar_label) ? defaults.colorbar_label :
                         colorbar_label,
        plotContour = plotContour,
        colormap = colormap,
        kwargs...)
end

function plot_glacier_dem(glacier::Glacier2D;
        title::Union{Nothing, String} = nothing,
        colorbar_label::Union{Nothing, String} = nothing,
        plotContour::Bool = true,
        colormap = :terrain,
        scale_text_size::Union{Nothing, Float64} = nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        logPlot::Bool = false)
    @assert !isempty(glacier.S) "The glacier DEM is empty."

    has_lonlat = haskey(glacier.Coords, "lon") && haskey(glacier.Coords, "lat") &&
                 !isempty(glacier.Coords["lon"]) && !isempty(glacier.Coords["lat"])
    x = has_lonlat ? glacier.Coords["lon"] : glacier.Δx .* collect(1.0:1.0:glacier.nx)
    y = has_lonlat ? glacier.Coords["lat"] : glacier.Δy .* collect(1.0:1.0:glacier.ny)

    dem_size = size(glacier.S)
    if !isempty(glacier.H₀) && size(glacier.H₀) == dem_size
        contour_mask = glacier.H₀ .> 0.0
    elseif !isempty(glacier.mask) && size(glacier.mask) == dem_size
        contour_mask = .!glacier.mask
    else
        contour_mask = trues(dem_size)
    end

    defaults = _dem_defaults(glacier)

    return _plot_gridded_data_core(glacier.S;
        lon = glacier.cenlon,
        lat = glacier.cenlat,
        x = x,
        y = y,
        rgi_id = glacier.rgi_id,
        Δx = glacier.Δx,
        mask = contour_mask,
        scale_text_size = scale_text_size,
        figsize = figsize,
        plotContour = plotContour,
        colormap = colormap,
        logPlot = logPlot,
        title = isnothing(title) ? defaults.title : title,
        colorbar_label = isnothing(colorbar_label) ? defaults.colorbar_label :
                         colorbar_label)
end

"""
    save_figure(fig, path)

Save `fig` to `path`, creating any missing parent directories automatically.
Returns `path`.
"""
function save_figure(fig::Figure, path::AbstractString)
    mkpath(dirname(path))
    CairoMakie.save(path, fig)
    return path
end
