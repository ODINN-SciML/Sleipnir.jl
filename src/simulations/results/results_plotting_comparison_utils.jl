###############################################################################
# Model-versus-observation plotting utilities (CairoMakie)
#
# Layout:
#   1. Internal helpers:  _finite_stats, _nearest_time_index, _observation_series,
#                         _observed_pixels
#   2. Public API:        plot_glacier_vs_observations, plot_field_histogram,
#                         plot_observations, plot_observation_timeline
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 1. Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _finite_stats(d) -> (; rmse, bias, n)

Root mean square error, mean bias and sample count over the finite entries of `d`.
Returns `NaN` statistics when nothing is finite, rather than throwing.
"""
function _finite_stats(d::AbstractMatrix)
    vals = filter(isfinite, vec(d))
    isempty(vals) && return (; rmse = NaN, bias = NaN, n = 0)
    return (; rmse = sqrt(sum(abs2, vals) / length(vals)),
        bias = sum(vals) / length(vals), n = length(vals))
end

"""
    _nearest_time_index(t, t_obs; atol) -> Int

Index of the simulated time closest to an observation time.

Warns, but does not throw, when the closest time is more than `atol` away: the survey can
still be plotted, but it is not matched in time.
"""
function _nearest_time_index(t::AbstractVector, t_obs::Real; atol::Real = 0.1)
    idx = argmin(abs.(t .- t_obs))
    gap = abs(t[idx] - t_obs)
    if gap > atol
        @warn "plot_glacier_vs_observations: no simulated state within $(atol) yr of the \
               observation" t_obs t_model=t[idx] gap
    end
    return idx
end

"""
    _observation_series(results, glacier, variable) -> (times, observations, label, unit)

Observation times and fields for `variable`, with the metadata used to label the panels.

`:H` comes from `glacier.thicknessData`, which has the survey dates (`results.H_ref` has no
times). `:V` comes from `results.V_ref` with `results.date_Vref`.
"""
function _observation_series(results::Results, glacier, variable::Symbol)
    if variable === :H
        td = glacier.thicknessData
        (isnothing(td) || isnothing(td.t) || isnothing(td.H)) && throw(ArgumentError(
            "plot_glacier_vs_observations(:H): glacier has no thicknessData with times."))
        return (td.t, td.H, "Ice thickness", "m")
    elseif variable === :V
        (isempty(results.V_ref) || isempty(first(results.V_ref))) && throw(ArgumentError(
            "plot_glacier_vs_observations(:V): results carry no reference velocities."))
        return (results.date_Vref, results.V_ref, "Surface velocity", "m yr⁻¹")
    else
        throw(ArgumentError("variable must be :H or :V, got $(variable)"))
    end
end

"""
    _observed_pixels(glacier, fields, times; max_points) -> (; x, y, v, t)

Position in km, value and time of the cells with an observation (zero or NaN means none).
North is up, as in the heatmaps. Only every k-th cell is kept, to have at most `max_points`.
"""
function _observed_pixels(glacier::Glacier2D, fields, times; max_points = 20_000)
    x, y = _glacier_xy(glacier)
    px = (; x = Float64[], y = Float64[], v = Float64[], t = Float64[])
    for (t, field) in zip(times, fields)
        if size(field) != (glacier.nx, glacier.ny)
            @warn "plot_observations: an observation at $(round(t; digits = 2)) is not on the grid of the glacier and is not plotted."
            continue
        end
        data = reverseForHeatmap(field, x, y)
        for c in findall(>(0), data)
            push!(px.x, (c[1] - 1) * glacier.Δx / 1000)
            push!(px.y, (c[2] - 1) * glacier.Δy / 1000)
            push!(px.v, data[c])
            push!(px.t, t)
        end
    end
    keep = 1:max(1, cld(length(px.v), max_points)):length(px.v)
    return map(a -> a[keep], px)
end

# ─────────────────────────────────────────────────────────────────────────────
# 2. Public API
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_glacier_vs_observations(
        results::Results,
        glacier,
        variable::Symbol;
        aggregate::Union{Nothing, Symbol} = nothing,
        mask::Union{Nothing, BitMatrix} = nothing,
        colormap = :viridis,
        diff_colormap = :redsblues,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        scale_text_size::Union{Nothing, Float64} = nothing,
        plotContour::Bool = false,
        title::Union{Nothing, String} = nothing,
    )

Compare a simulated field with its observations: modelled, observed and difference.

One row of three panels per observation time. The modelled panel is the simulated state at
the saved time closest to the observation, not the last time step.

# Arguments

  - `results::Results`: Simulation output, with the modelled field and the grid metadata.
  - `glacier`: The glacier, needed for `:H` because the survey dates are in `thicknessData`.
  - `variable::Symbol`: `:H` (ice thickness) or `:V` (surface velocity magnitude).
  - `aggregate::Union{Nothing,Symbol}`: `nothing` (default) plots every observation time.
    `:mean` gives a single row: for each cell, the observations and the modelled field are
    averaged over the times where this cell is observed.
  - `mask::Union{Nothing,BitMatrix}`: Cells to keep. Defaults to all the cells. In any case,
    only the cells with an observation are compared: zero means that there is no observation,
    as in `ThicknessData` and `SurfaceVelocityData`.
  - `colormap`, `diff_colormap`, `figsize`, `scale_text_size`, `plotContour`, `title`:
    Optional plotting parameters.

# Returns

  - A `CairoMakie.Figure`.

# Notes

  - Modelled and observed panels share a colour range so the two are directly comparable;
    the difference (model minus observation) panel is symmetric about zero.
  - Each difference panel is titled with its RMSE and mean bias over the masked cells.
  - The observations more than 0.1 yr outside the simulated period are not plotted.
"""
function plot_glacier_vs_observations(
        results::Results,
        glacier,
        variable::Symbol;
        aggregate::Union{Nothing, Symbol} = nothing,
        mask::Union{Nothing, BitMatrix} = nothing,
        colormap = :viridis,
        diff_colormap = :redsblues,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        scale_text_size::Union{Nothing, Float64} = nothing,
        plotContour::Bool = false,
        title::Union{Nothing, String} = nothing
)
    (isnothing(aggregate) || aggregate === :mean) ||
        throw(ArgumentError("aggregate must be nothing or :mean, got $(aggregate)."))
    meta = _results_plot_metadata(results; plotContour = plotContour)
    t_obs, obs, label, unit = _observation_series(results, glacier, variable)
    modelled_series = variable === :H ? results.H : results.V
    (isempty(modelled_series) || isempty(first(modelled_series))) && throw(ArgumentError(
        "plot_glacier_vs_observations(:$(variable)): results carry no simulated $(variable)."))
    length(t_obs) == length(obs) || throw(ArgumentError(
        "plot_glacier_vs_observations(:$(variable)): the observations have no dates."))

    isempty(results.t) && throw(ArgumentError(
        "plot_glacier_vs_observations(:$(variable)): the results have no times."))
    # Only the observations inside the simulated period are compared
    inside = findall(tᵢ -> results.t[begin] - 0.1 <= tᵢ <= results.t[end] + 0.1, t_obs)
    length(inside) < length(t_obs) &&
        @info "plot_glacier_vs_observations: $(length(t_obs) - length(inside)) observations outside the simulated period are not plotted."
    isempty(inside) && throw(ArgumentError(
        "plot_glacier_vs_observations(:$(variable)): no observation in the simulated period."))
    t_obs, obs = t_obs[inside], obs[inside]

    # Pair each observation with the simulated state at the nearest saved time
    idxs = [_nearest_time_index(results.t, tᵢ) for tᵢ in t_obs]
    rows = if isnothing(aggregate)
        [(t_obs[i], modelled_series[idxs[i]], obs[i]) for i in eachindex(t_obs)]
    else
        # Average both sides cell by cell, over the times where the cell is observed
        seen = [o .> 0 for o in obs]
        n = max.(sum(seen), 1)
        obs_mean = sum(ifelse.(s, o, 0.0) for (o, s) in zip(obs, seen)) ./ n
        mod_mean = sum(ifelse.(s, modelled_series[i], 0.0) for (i, s) in zip(idxs, seen)) ./
                   n
        [(sum(t_obs) / length(t_obs), mod_mean, obs_mean)]
    end

    base_mask = isnothing(mask) ? trues(size(first(rows)[3])) : mask
    panel_size = 360 # pixels of the longest side of a map
    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)
    fig = Figure(; figKwargs...)

    for (row, (tᵢ, mod_raw, obs_raw)) in enumerate(rows)
        # Zero (and NaN) means no observation, so the mask is per row
        cells = base_mask .& (obs_raw .> 0)
        modelled = float.(mod_raw)
        observed = float.(obs_raw)
        modelled[.!cells] .= NaN
        observed[.!cells] .= NaN
        difference = modelled .- observed
        stats = _finite_stats(difference)

        shared = filter(isfinite, vcat(vec(modelled), vec(observed)))
        if isempty(shared)
            @warn "plot_glacier_vs_observations: no observed cell left in the mask at $(round(tᵢ; digits = 2)), this row is not plotted."
            continue
        end
        lo, hi = extrema(shared)
        hi > lo || (hi = lo + 1)
        dmax = maximum(abs, filter(isfinite, vec(difference)); init = 0.0)
        dmax = dmax > 0 ? dmax : 1.0

        when = isnothing(aggregate) ? "$(round(tᵢ; digits = 2))" :
               "mean over $(length(idxs)) obs"
        panels = (
            ("Modelled ($when)", modelled, colormap, (lo, hi)),
            ("Observed ($when)", observed, colormap, (lo, hi)),
            (
                @sprintf("Difference\nRMSE %.3g, bias %+.3g %s", stats.rmse, stats.bias,
                    unit),
                difference, diff_colormap, (-dmax, dmax))
        )

        for (col, (paneltitle, data, cmap, crange)) in enumerate(panels)
            nx, ny = size(data)
            # Size of the panel from the shape of the grid, so that the map is large with
            # respect to the labels. The cells are not interpolated: the grid can be coarse.
            w, h = panel_size .* (nx, ny) ./ max(nx, ny)
            ax = Axis(fig[row, 2col - 1], aspect = DataAspect(), width = w, height = h,
                title = paneltitle, titlesize = 13)
            hm = heatmap!(ax, reverseForHeatmap(data, results.x, results.y),
                colormap = cmap, colorrange = crange, interpolate = false)
            Colorbar(fig[row, 2col], hm; label = col == 3 ? "Δ $unit" : unit, height = h)
            plotContour && _overlay_contour!(ax, meta.contour)
            _decorate_geo_axis!(ax, nx, ny, meta.lon, meta.lat, meta.Δx;
                scale_text_size = something(scale_text_size, 12.0), num_vars = 3)
        end
    end

    fig_title = isnothing(title) ? "$label — $(meta.rgi_id)" : "$title — $(meta.rgi_id)"
    fig[0, :] = Label(fig, fig_title, fontsize = 14, font = :bold)
    isnothing(figsize) && resize_to_layout!(fig)
    return fig
end

"""
    plot_field_histogram(
        data::AbstractMatrix;
        references = (),
        mask::Union{Nothing, BitMatrix} = nothing,
        bins::Int = 50,
        logScale::Bool = false,
        xlabel::String = "value",
        title::Union{Nothing, String} = nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        color = (:steelblue, 0.7),
    )

Histogram of an inverted field, with reference values (for example from the literature)
drawn as vertical lines, to see if the inverted values are in a plausible range.

# Arguments

  - `data::AbstractMatrix`: The field. Non-finite entries are dropped.
  - `references`: Pairs of `label => value` drawn as labelled vertical lines, e.g.
    `["Cuffey & Paterson 0°C" => 7.6e-17, "−10°C" => 1.1e-17]`.
  - `mask::Union{Nothing,BitMatrix}`: Cells to include. Use it for masked fields: the cells
    without ice keep their initial value and make a false peak.
  - `bins::Int`: Number of histogram bins.
  - `logScale::Bool`: Bin in `log10` space, for fields spanning orders of magnitude.
  - `xlabel`, `title`, `figsize`, `color`: Optional plotting parameters.

# Returns

  - A `CairoMakie.Figure`.
"""
function plot_field_histogram(
        data::AbstractMatrix;
        references = (),
        mask::Union{Nothing, BitMatrix} = nothing,
        bins::Int = 50,
        logScale::Bool = false,
        xlabel::String = "value",
        title::Union{Nothing, String} = nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        color = (:steelblue, 0.7)
)
    selected = isnothing(mask) ? vec(data) : vec(data[mask])
    vals = collect(filter(isfinite, selected))
    if logScale
        vals = filter(v -> v > 0, vals)
    end
    @assert !isempty(vals) "plot_field_histogram: no finite values left after masking."

    xs = logScale ? log10.(vals) : vals
    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}(:size => (700, 420)) :
                Dict{Symbol, Any}(:size => figsize)
    fig = Figure(; figKwargs...)
    ax = Axis(fig[1, 1],
        xlabel = logScale ? "log₁₀($xlabel)" : xlabel,
        ylabel = "count",
        title = isnothing(title) ? "" : title)

    hist!(ax, xs, bins = bins, color = color)

    # Reference lines, cycling colours so several stay distinguishable
    palette = [:firebrick, :darkorange, :seagreen, :purple, :black]
    for (i, ref) in enumerate(references)
        reflabel, refvalue = ref isa Pair ? (first(ref), last(ref)) : ("", ref)
        (logScale && refvalue <= 0) && continue
        pos = logScale ? log10(refvalue) : refvalue
        vlines!(ax, [pos], color = palette[mod1(i, length(palette))],
            linestyle = :dash, linewidth = 2, label = String(reflabel))
    end
    isempty(references) || axislegend(ax, position = :rt, framevisible = false,
        labelsize = 10)

    # Median of the field itself, so the comparison with the references is quantitative
    med = Statistics.median(xs)
    vlines!(ax, [med], color = :grey20, linewidth = 2)
    text!(ax, @sprintf(" median %.3g", logScale ? 10^med : med),
        position = (med, 0.0), align = (:left, :bottom), fontsize = 10, color = :grey20)

    resize_to_layout!(fig)
    return fig
end

"""
    plot_observations(glacier::Glacier2D; tspan = nothing, max_points = 20_000)

Make a 3D figure of where and when a glacier has observations: x and y are in km (north up)
and z is the time in decimal years. Ice thickness surveys and surface velocity maps are
drawn at their dates, coloured by value. Zero means no observation, so those cells are not
drawn. The glacier-wide `:dhdt` is only in [`plot_observation_timeline`](@ref).

  - `tspan`: Time window, drawn as two red rectangles, for example `glacier_tspan(params, glacier)`.
  - `max_points`: Maximum pixels drawn for each kind. If there are more, every k-th is drawn.
"""
function plot_observations(
        glacier::Glacier2D; tspan = nothing, max_points::Int = 20_000)
    td, vd = glacier.thicknessData, glacier.velocityData
    H = _observed_pixels(glacier, isnothing(td) ? [] : something(td.H, []), tdata(td);
        max_points)
    V = _observed_pixels(glacier, isnothing(vd) ? [] : something(vd.vabs, []),
        tdata(vd, MeanDateVelocityMapping()); max_points)
    (isempty(H.v) && isempty(V.v)) && throw(ArgumentError(
        "plot_observations: glacier $(glacier.rgi_id) has no observed pixels."))

    x_km = (glacier.nx - 1) * glacier.Δx / 1000
    y_km = (glacier.ny - 1) * glacier.Δy / 1000
    fig = Figure(size = (900, 650))
    ax = Axis3(fig[1, 1]; xlabel = "x (km)", ylabel = "y (km)", zlabel = "Time (year)",
        title = "Observations — $(glacier.rgi_id)", aspect = (x_km, y_km, max(x_km, y_km)),
        azimuth = 1.4π, elevation = π / 8)
    for (px, colormap, label) in (
        (H, :YlGnBu, "Ice thickness (m)"), (V, :viridis, "Surface velocity (m yr⁻¹)"))
        isempty(px.v) && continue
        sc = scatter!(ax, px.x, px.y, px.t; color = px.v, colormap, markersize = 6)
        Colorbar(fig[1, end + 1], sc; label)
    end

    # Rectangles at the start and the end of the time window
    foreach(
        z -> lines!(ax, Point3f.([0, x_km, x_km, 0, 0], [0, 0, y_km, y_km, 0], z);
            color = :firebrick),
        something(tspan, ()))
    return fig
end

"""
    plot_observation_timeline(glaciers::Vector{<:Glacier2D}; tspan = nothing, glacier_tspans = Dict())

Plot when each glacier has observations, one row per glacier. The marks come from
[`observations`](@ref): a dot at the date of each observation (`:H` blue, `:V` orange,
`:dhdt` red) and a line for its time span, if it has one. `glacier_tspans` is a grey bar for
each glacier that has an entry by RGI ID, for example `params.simulation.glacier_tspans`.
`tspan` is drawn as two vertical dashed lines. Returns a `CairoMakie.Figure`.
"""
function plot_observation_timeline(
        glaciers::Vector{<:Glacier2D}; tspan = nothing, glacier_tspans = Dict())
    @assert !isempty(glaciers) "plot_observation_timeline: no glaciers to plot."
    n = length(glaciers)
    colors = Dict(:H => :steelblue, :V => :darkorange, :dhdt => :firebrick)
    obs = [(; o..., row = i) for (i, g) in enumerate(glaciers) for o in observations(g)]
    spans = filter(o -> o.t_end > o.t_start, obs)
    windows = [Rect(t[1], i - 0.4, t[2] - t[1], 0.8)
               for (i, g) in enumerate(glaciers) if haskey(glacier_tspans, g.rgi_id)
               for t in (glacier_tspans[g.rgi_id],)]

    fig = Figure(size = (900, 150 + 35n))
    ax = Axis(fig[1, 1]; title = "Observations (H blue, V orange, dhdt red)",
        xlabel = "Time (year)", yticks = (1:n, [g.rgi_id for g in glaciers]),
        yreversed = true, limits = (nothing, nothing, 0.4, n + 0.6))
    isempty(windows) || poly!(ax, windows; color = (:grey, 0.25))
    isempty(spans) || rangebars!(ax, [o.row for o in spans], [o.t_start for o in spans],
        [o.t_end for o in spans]; direction = :x, linewidth = 3,
        color = [colors[o.kind] for o in spans])
    scatter!(ax, [o.t for o in obs], [o.row for o in obs];
        color = [colors[o.kind] for o in obs], markersize = 9)
    isnothing(tspan) ||
        vlines!(ax, collect(Float64, tspan); color = :black, linestyle = :dash)
    return fig
end
