###############################################################################
# Model-versus-observation plotting utilities (CairoMakie)
#
# Layout:
#   1. Internal helpers:  _finite_stats, _nearest_time_index, _observation_series
#   2. Public API:        plot_glacier_vs_observations, plot_field_histogram
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 1. Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _finite_stats(diff) -> (; rmse, bias, n)

Root mean square error, mean bias and sample count over the finite entries of `diff`.
Returns `NaN` statistics when nothing is finite, rather than throwing.
"""
function _finite_stats(diff::AbstractMatrix)
    vals = filter(isfinite, vec(diff))
    isempty(vals) && return (; rmse = NaN, bias = NaN, n = 0)
    return (; rmse = sqrt(sum(abs2, vals) / length(vals)),
        bias = sum(vals) / length(vals), n = length(vals))
end

"""
    _nearest_time_index(t, t_obs; atol) -> Int

Index of the simulated time closest to an observation time.

Emits a warning past `atol` instead of throwing: a survey slightly outside the saved times
is still worth looking at, it just must not be silently presented as time matched.
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

`:H` comes from `glacier.thicknessData`, the only place carrying survey dates — `Results`
stores `H_ref` with no time vector, so plotting against it cannot be time matched. `:V`
comes from `results.V_ref` with `results.date_Vref`.
"""
function _observation_series(results::Results, glacier, variable::Symbol)
    if variable === :H
        td = glacier.thicknessData
        (isnothing(td) || isnothing(td.t) || isnothing(td.H)) && throw(ArgumentError(
            "plot_glacier_vs_observations(:H): glacier has no thicknessData with times."))
        return (td.t, td.H, "Ice thickness", "m")
    elseif variable === :V
        isempty(results.V_ref) && throw(ArgumentError(
            "plot_glacier_vs_observations(:V): results carry no reference velocities."))
        return (results.date_Vref, results.V_ref, "Surface velocity", "m yr⁻¹")
    else
        throw(ArgumentError("variable must be :H or :V, got $(variable)"))
    end
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

Compare a simulated field against its observations: modelled, observed and difference.

One row of three panels per observation time. The modelled panel is the simulated state at
the saved time nearest the observation, **not** the last time step — comparing an inverted
field against an observation from a different date is the most common way to misread an
inversion, so the pairing is done here rather than left to the caller.

# Arguments

  - `results::Results`: Simulation output, supplying the modelled field and the grid metadata.
  - `glacier`: The glacier, needed for `:H` because survey dates live in `thicknessData`.
  - `variable::Symbol`: `:H` (ice thickness) or `:V` (surface velocity magnitude).
  - `aggregate::Union{Nothing,Symbol}`: `nothing` (default) plots every observation time.
    `:mean` averages the observations and the modelled field over the observation window
    into a single row, matching what a window-averaged velocity loss actually compares.
  - `mask::Union{Nothing,BitMatrix}`: Cells to keep. Defaults to all the cells. In any case,
    only the cells with an observation are compared: zero means that there is no observation,
    as in `ThicknessData` and `SurfaceVelocityData`.
  - `colormap`, `diff_colormap`, `figsize`, `scale_text_size`, `plotContour`, `title`:
    Optional plotting parameters.

# Returns

  - A `CairoMakie.Figure`.

# Notes

  - Modelled and observed panels share a colour range so the two are directly comparable;
    the difference panel is symmetric about zero.
  - Each difference panel is titled with its RMSE and mean bias over the masked cells.
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
    @assert isnothing(aggregate) || aggregate === :mean "aggregate must be nothing or :mean"
    meta = _results_plot_metadata(results; plotContour = plotContour)
    t_obs, obs, label, unit = _observation_series(results, glacier, variable)
    modelled_series = variable === :H ? results.H : results.V
    @assert length(t_obs)==length(obs) "observation times and fields disagree in length"
    @assert !isempty(t_obs) "no observations to compare against"

    # Pair each observation with the simulated state at the nearest saved time
    idxs = [_nearest_time_index(results.t, tᵢ) for tᵢ in t_obs]
    pairs = if isnothing(aggregate)
        [(t_obs[i], modelled_series[idxs[i]], obs[i]) for i in eachindex(t_obs)]
    else
        # Average both sides over the same set of observation times
        mod_mean = sum(modelled_series[i] for i in idxs) / length(idxs)
        obs_mean = sum(obs) / length(obs)
        [(sum(t_obs) / length(t_obs), mod_mean, obs_mean)]
    end

    base_mask = isnothing(mask) ? trues(size(first(pairs)[3])) : mask
    panel_size = 360 # pixels of the longest side of a map
    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)
    fig = Figure(; figKwargs...)

    for (row, (tᵢ, mod_raw, obs_raw)) in enumerate(pairs)
        # Zero (and NaN) means no observation, so the mask is per row
        cells = base_mask .& (obs_raw .> 0)
        modelled = copy(float.(mod_raw))
        observed = copy(float.(obs_raw))
        modelled[.!cells] .= NaN
        observed[.!cells] .= NaN
        difference = modelled .- observed
        stats = _finite_stats(difference)

        shared = filter(isfinite, vcat(vec(modelled), vec(observed)))
        @assert !isempty(shared) "no finite cells left where model and observation overlap"
        lo, hi = extrema(shared)
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
    resize_to_layout!(fig)
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

Distribution of an inverted field with labelled reference values marked on it.

The point is to read a recovered parameter against the values physics expects: a histogram
of `A` or `C` against literature values tells you at a glance whether the inversion landed
in a plausible range or merely converged. Reference values are supplied by the caller, so
this stays a rendering function with no physics in it.

# Arguments

  - `data::AbstractMatrix`: The field. Non-finite entries are dropped.
  - `references`: Pairs of `label => value` drawn as labelled vertical lines, e.g.
    `["Cuffey & Paterson 0°C" => 7.6e-17, "−10°C" => 1.1e-17]`.
  - `mask::Union{Nothing,BitMatrix}`: Cells to include. **Pass this for any masked field.**
    Ice-free cells keep their seed value forever and pile up as a spurious spike.
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
