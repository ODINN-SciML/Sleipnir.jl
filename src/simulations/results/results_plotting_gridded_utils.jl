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
        logPlot = false,
        title::Union{Nothing, String} = nothing,
        colorbar_label::Union{Nothing, String} = nothing
) where {F <: AbstractFloat}
    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)

    nx, ny = size(mask)
    ctr = plotContour ? Contour.contour(collect(1:nx), 1+ny .- collect(1:ny), mask, 0.5) :
          nothing

    @assert !isempty(gridded_data) "There is no data."
    if typeof(gridded_data) <: Vector
        @assert length(gridded_data)>0 "Data is an empty vector"
        @assert (isnothing(timeIdx)) || (size(gridded_data, 1)>=timeIdx) "The provided index=$(timeIdx) is greater than the size of the vector which is $(size(gridded_data,1))"
        current_data = isnothing(timeIdx) ? gridded_data[end] : gridded_data[timeIdx]
    else
        current_data = gridded_data
    end

    @assert size(current_data) == size(mask) "The data and mask must have the same size."

    current_data_masked = copy(current_data)
    current_data_masked[.!mask] .= NaN

    if logPlot
        global_min = lower_bound(current_data_masked)
        global_max = maximum(replace(current_data_masked, NaN => 0.0))
    else
        global_min, global_max = finite_extrema(current_data_masked)
    end

    figKwargs[:layout] = GridLayout(2, 2)
    fig = Figure(; figKwargs...)

    ax_row = 1
    ax_col = 1
    ax = Axis(fig[ax_row, ax_col], aspect = DataAspect())
    data = deepcopy(gridded_data)

    if typeof(data) <: Vector
        @assert length(data)>0 "Data is an empty vector"
        @assert (isnothing(timeIdx)) || (size(data, 1)>=timeIdx) "The provided index=$(timeIdx) is greater than the size of the vector which is $(size(data,1))"
        data = isnothing(timeIdx) ? data[end] : data[timeIdx]
    end

    nx, ny = size(data)

    data[.!mask] .= NaN

    hm = heatmap!(ax, reverseForHeatmap(data, x, y), colormap = colormap,
        colorrange = (logPlot ? global_min*0.85 : global_min, global_max),
        colorscale = logPlot ? log10 : identity)
    cb_kwargs = isnothing(colorbar_label) ? NamedTuple() : (; label = colorbar_label)
    cb = Colorbar(fig[ax_row, ax_col + 1], hm; cb_kwargs...)
    Observables.connect!(cb.height, @lift CairoMakie.Fixed($(viewport(ax.scene)).widths[2]))

    if plotContour
        for curve in ctr.lines
            xs = first.(curve.vertices)
            ys = last.(curve.vertices)
            lines!(ax, xs, ys, color = :black, linewidth = 1)
        end
    end

    ax.xlabel = "Longitude (°)"
    ax.ylabel = "Latitude (°)"
    ax.xticks=([round(nx/2)], ["$(round(lon;digits=6))"])
    ax.yticks=([round(ny/2)], ["$(round(lat;digits=6))"])
    ax.yticklabelrotation = π/2
    ax.ylabelpadding = 5
    ax.yticklabelalign = (:center, :bottom)

    scale_width = 0.10*nx
    scale_number = round(Δx * scale_width / 1000; digits = 1) # Convert to km
    textsize = isnothing(scale_text_size) ? 1.2*scale_width : scale_text_size

    poly!(ax, Rect(nx - round(0.15*nx), round(0.075*ny), scale_width, scale_width/10),
        color = :black)
    text!(ax,
        "$scale_number km",
        position = (nx - round(0.15*nx) + scale_width/16, round(0.075*ny) + scale_width/10),
        fontsize = textsize)

    fig_title = isnothing(title) ? rgi_id : "$title — $rgi_id"
    fig[0, :] = Label(fig, fig_title, fontsize = 14, font = :bold)
    resize_to_layout!(fig)
    return fig
end

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

Plot cumulative mass-balance map from MB fields stored at each forward MB callback.

# Arguments

    - `results::Results`: Results containing callback MB maps in `results.MB`.
    - `kwargs...`: Keyword arguments forwarded to `plot_cumulative_gridded_data`.

# Returns

    - `Figure`: Cumulative MB figure.
"""
function plot_cumulative_mb(results::Results;
        title::String = "Cumulative Mass Balance",
        colorbar_label::String = "m w.e.",
        kwargs...)
    if isempty(results.MB) || isempty(results.MB[begin])
        @warn "No mass balance callback history in results; skipping cumulative MB plot. " *
              "Make sure the simulation was run with use_MB=true."
        return nothing
    end
    return plot_cumulative_gridded_data(results.MB, results;
        title = title, colorbar_label = colorbar_label, kwargs...)
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
