###############################################################################
# Plots to inspect a simulation or an inversion against the observations (CairoMakie)
#
# Layout:
#   1. Style:          paper_theme
#   2. Maps:           MapPanel, plot_field_maps, plot_surface_change, plot_bed_maps
#   3. Observations:   plot_fit_to_observations
#
# Rules of all the figures: the colour range is the full range of the data (nothing is clipped),
# differences are symmetric about zero, and the cells without data are light grey, never white.
# They only need a `Results` and a `Glacier2D`. Save them as PDF for a paper (`save("fig.pdf", fig)`).
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 1. Style
# ─────────────────────────────────────────────────────────────────────────────

"""
    paper_theme() -> Theme

Theme for the figures of papers and reports: small sans-serif text, thin frames, no grid and
thin colorbars. The plot functions of this file use it already. For your own plots, use
`with_theme(paper_theme()) do ... end`.
"""
function paper_theme()
    return Theme(
        fontsize = 10,
        Axis = (titlesize = 10, titlefont = :regular, titlealign = :left,
            xlabelsize = 10, ylabelsize = 10, xticklabelsize = 9, yticklabelsize = 9,
            spinewidth = 0.6, xtickwidth = 0.6, ytickwidth = 0.6,
            xgridvisible = false, ygridvisible = false),
        Colorbar = (labelsize = 9, ticklabelsize = 8, width = 9, spinewidth = 0.5,
            tickwidth = 0.5),
        Legend = (framevisible = false, labelsize = 9, patchsize = (14, 8)),
        Lines = (linewidth = 1.3,)
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# 2. Maps
# ─────────────────────────────────────────────────────────────────────────────

"""
    MapPanel(data, title; cmap = :viridis, label = "", symmetric = false, limits = nothing, nodata = nothing)

One map of [`plot_field_maps`](@ref).

  - `data`: the field on the grid of the glacier. `NaN` cells are not drawn.
  - `title`: title of the map.
  - `cmap`: colormap.
  - `label`: label of the colorbar, for example `"thickness (m)"`.
  - `symmetric`: if `true`, the colour range is symmetric about zero (for differences).
  - `limits`: colour range `(min, max)`. By default it is the full range of `data`.
  - `nodata`: cells that have no data, drawn in light grey. White is the colour of zero.
"""
Base.@kwdef struct MapPanel
    data::Matrix{Float64}
    title::String
    cmap::Any = :viridis
    label::String = ""
    symmetric::Bool = false
    limits::Union{Nothing, Tuple{Float64, Float64}} = nothing
    nodata::Union{Nothing, BitMatrix} = nothing
end
function MapPanel(data::AbstractMatrix, title::AbstractString; kwargs...)
    return MapPanel(; data = Float64.(data), title = String(title), kwargs...)
end

# The map in the orientation of the heatmap
function _orient(A, results::Results)
    return reverseForHeatmap(Matrix{eltype(results.x)}(A), results.x, results.y)
end

# A round length in km for the scale bar, about a fifth of the width of the map
function _scalebar_length_km(width_km)
    options = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
    return options[argmin(abs.(options .- width_km / 5))]
end

function _scalebar!(ax, nx, ny, Δx)
    km = _scalebar_length_km(nx * Δx / 1000)
    len = km * 1000 / Δx
    x0, y0 = nx - 0.04 * nx - len, 0.05 * ny
    lines!(ax, [x0, x0 + len], [y0, y0]; color = :black, linewidth = 2.5)
    text!(ax, x0 + len / 2, y0 + 0.015 * ny; text = "$(round(km; digits = 1)) km",
        align = (:center, :bottom), fontsize = 8)
end

# Colour range of a panel: the full range of the data
function _panel_limits(panel::MapPanel)
    isnothing(panel.limits) || return panel.limits
    v = filter(isfinite, panel.data)
    isempty(v) && return (0.0, 1.0)
    if panel.symmetric
        m = maximum(abs, v)
        return m > 0 ? (-m, m) : (-1.0, 1.0)
    end
    lo, hi = extrema(v)
    return hi > lo ? (lo, hi) : (lo, lo + 1.0)
end

function _map_panel!(pos, panel::MapPanel, results::Results, outline)
    nx, ny = size(panel.data)
    ax = Axis(pos[1, 1]; title = panel.title, aspect = DataAspect())
    hidedecorations!(ax)
    isnothing(panel.nodata) ||
        heatmap!(ax, _orient(ifelse.(panel.nodata, 1.0, NaN), results);
            colormap = [:lightgray, :lightgray], colorrange = (0, 1), interpolate = false)
    hm = heatmap!(ax, _orient(panel.data, results); colormap = panel.cmap,
        colorrange = _panel_limits(panel), interpolate = false)
    Colorbar(pos[1, 2], hm; label = panel.label)
    isnothing(outline) ||
        linesegments!(ax, _mask_outline(_orient(float.(outline), results) .> 0);
            color = :black, linewidth = 0.6)
    _scalebar!(ax, nx, ny, results.Δx)
    return ax
end

"""
    plot_field_maps(results, panels; ncols = 3, outline = results.H[begin] .> 0, title = nothing) -> Figure

Maps of fields of a glacier in one figure, with the same style: one colorbar for each map, a
scale bar, the outline of the glacier, and the letters (a), (b), ... in the titles. This is the
function that the other map figures of this file are made with; use it for any other field.

  - `panels`: a vector of [`MapPanel`](@ref).
  - `ncols`: number of columns.
  - `outline`: cells inside the black outline. By default the cells with ice at the start.
  - `title`: title of the figure.
"""
function plot_field_maps(
        results::Results,
        panels::AbstractVector{MapPanel};
        ncols::Integer = 3,
        outline::Union{Nothing, AbstractMatrix{Bool}} = results.H[begin] .> 0,
        title::Union{Nothing, AbstractString} = nothing
)
    isempty(panels) && throw(ArgumentError("plot_field_maps needs at least one MapPanel."))
    return with_theme(paper_theme()) do
        n = length(panels)
        ncols = min(ncols, n)
        nrows = cld(n, ncols)
        nx, ny = size(first(panels).data)
        height_map = 280 * ny / nx + 30
        fig = Figure(size = (
            ncols * 330, nrows * height_map + (isnothing(title) ? 10 : 40)))
        isnothing(title) || Label(fig[0, 1:ncols], title; fontsize = 12, font = :bold)
        for (k, panel) in enumerate(panels)
            row, col = fld1(k, ncols), mod1(k, ncols)
            lettered = MapPanel(;
                data = panel.data, title = "($('a' + k - 1)) $(panel.title)",
                cmap = panel.cmap, label = panel.label, symmetric = panel.symmetric,
                limits = panel.limits, nodata = panel.nodata)
            _map_panel!(fig[row, col], lettered, results, outline)
        end
        colgap!(fig.layout, 8)
        rowgap!(fig.layout, 8)
        fig
    end
end

# Cells where ice is allowed (the outline and its buffer), or all the cells
function _free_cells(glacier, A)
    mask = glacier.mask
    return size(mask) == size(A) ? .!mask : trues(size(A))
end

# Surface of the model at time t
function _model_surface(results::Results, B, t)
    return B .+ max.(results.H[argmin(abs.(results.t .- t))], 0)
end

"""
    plot_surface_change(results, glacier; B = results.B, title = nothing) -> Figure

Change of the surface between the first and the last dated DEM of the glacier
(`glacier.surfaceData`) that are inside the simulated period: (a) from the DEMs, (b) from the
model, with the same colour range, and (c) the difference model minus DEMs.

Red is a lower surface (negative) and blue a higher one (positive). Cells where one of the two
DEMs has no usable value are light grey, for example the cells of the Copernicus DEM that were
filled from other DEMs and have no error and no date. The model is only compared where both
DEMs exist.

`B` is the bed of the model (by default `results.B`): pass the inverted bed when the bed was
inverted.
"""
function plot_surface_change(
        results::Results, glacier;
        B::AbstractMatrix = results.B,
        title::Union{Nothing, AbstractString} = nothing
)
    sd = glacier.surfaceData
    isnothing(sd) && throw(ArgumentError(
        "plot_surface_change: the glacier has no dated DEMs (surfaceData)."))
    inside = findall(t -> results.t[begin] <= t <= results.t[end], sd.t)
    length(inside) >= 2 || throw(ArgumentError(
        "plot_surface_change needs two dated DEMs inside the simulated period."))
    ka, kb = first(inside), last(inside)
    ta, tb = sd.t[ka], sd.t[kb]
    free = _free_cells(glacier, sd.S[ka])
    ok = isfinite.(sd.S[ka]) .& isfinite.(sd.S[kb]) .& free
    nodata = BitMatrix(free .& .!ok)
    dem = ifelse.(ok, sd.S[kb] .- sd.S[ka], NaN)
    model = ifelse.(ok, _model_surface(results, B, tb) .- _model_surface(results, B, ta), NaN)
    m = maximum(abs, filter(isfinite, vcat(vec(dem), vec(model))))
    limits = (-m, m)
    # Red is negative and blue positive
    cmap = Reverse(:balance)
    when = "$(_date_label(ta)) → $(_date_label(tb))"
    panels = [
        MapPanel(dem, "DEMs, $(when)"; cmap, label = "surface change (m)", limits, nodata),
        MapPanel(
            model, "Model, $(when)"; cmap, label = "surface change (m)", limits, nodata),
        MapPanel(model .- dem, "Model − DEMs"; cmap, label = "difference (m)",
            symmetric = true, nodata)
    ]
    return plot_field_maps(results, panels; ncols = 3,
        title = something(title, "Surface change of $(results.rgi_id)"))
end

# Local mean of A over the valid cells of a window of half-width r (NaN where no cell is valid)
function _local_mean(A, valid, r::Integer)
    out = fill(NaN, size(A))
    for j in axes(A, 2), i in axes(A, 1)

        valid[i, j] || continue
        s, n = 0.0, 0
        for dj in (-r):r, di in (-r):r

            ii, jj = i + di, j + dj
            (1 <= ii <= size(A, 1) && 1 <= jj <= size(A, 2) && valid[ii, jj]) || continue
            s += A[ii, jj]
            n += 1
        end
        out[i, j] = s / n
    end
    return out
end

# Cells of the outline whose 8 neighbours are also inside it. At the margin, the mean over a window
# is one-sided and the steep valley walls give a relief that is the slope, not a feature of the bed.
function _inner_cells(outline)
    inner = copy(BitMatrix(outline))
    for I in CartesianIndices(outline), dj in -1:1, di in -1:1
        J = I + CartesianIndex(di, dj)
        (checkbounds(Bool, outline, J) && outline[J]) || (inner[I] = false)
    end
    return inner
end

# Relief of a field at the scale of one cell: the field minus its mean over 3 × 3 cells, where the
# whole window is inside the outline (NaN at the margin)
function _cell_relief(A, outline)
    return ifelse.(_inner_cells(outline), A .- _local_mean(A, outline, 1), NaN)
end

# Slope (m/m) of a field with central differences
function _gradient_norm(A, Δx, Δy)
    gx = zeros(size(A))
    gy = zeros(size(A))
    gx[2:(end - 1), :] .= (A[3:end, :] .- A[1:(end - 2), :]) ./ (2Δx)
    gy[:, 2:(end - 1)] .= (A[:, 3:end] .- A[:, 1:(end - 2)]) ./ (2Δy)
    return hypot.(gx, gy)
end

"""
    plot_bed_maps(results, glacier; B = results.B, title = nothing) -> Figure

The bed and what the ice flow model needs from it: (a) the elevation of the bed, (b) the ice
thickness at the start, (c) the slope of the bed, and (d) the relief of the bed at the scale of
one cell (the bed minus its mean over 3 × 3 cells, not drawn at the margin of the glacier, where
the steep slope of the valley walls would give a relief). The ice flow model (SIA) assumes that the
thickness changes slowly over a few ice thicknesses, so a large relief at the scale of a cell is
a warning. `B` is the bed of the model, as in [`plot_surface_change`](@ref).
"""
function plot_bed_maps(
        results::Results, glacier;
        B::AbstractMatrix = results.B,
        title::Union{Nothing, AbstractString} = nothing
)
    free = _free_cells(glacier, B)
    outline = results.H[begin] .> 0
    masked(A, cells) = ifelse.(cells, A, NaN)
    H₀ = results.H[begin]
    relief = _cell_relief(B, outline)
    panels = [
        MapPanel(masked(B, free), "Bed elevation"; cmap = :terrain, label = "m a.s.l."),
        MapPanel(masked(H₀, outline), "Ice thickness at the start"; label = "m"),
        MapPanel(
            masked(_gradient_norm(B, results.Δx, results.Δy), free), "Slope of the bed";
            cmap = :magma, label = "m/m"),
        MapPanel(masked(relief, outline), "Bed relief at the cell scale (3 × 3)";
            cmap = :balance, label = "m", symmetric = true)
    ]
    return plot_field_maps(results, panels; ncols = 2,
        title = something(title, "Bed of $(results.rgi_id)"))
end

# ─────────────────────────────────────────────────────────────────────────────
# 3. Observations
# ─────────────────────────────────────────────────────────────────────────────

_rmse_bias(d) = (sqrt(mean(abs2, d)), mean(d))

function _scatter_fit!(
        pos, observed, modelled, title; xlabel, ylabel, color = (:black, 0.45),
        colormap = :viridis, colorrange = nothing)
    ax = Axis(pos; title, xlabel, ylabel, aspect = 1)
    kw = isnothing(colorrange) ? (;) : (; colorrange)
    sc = scatter!(ax, observed, modelled; markersize = 3.5, color, colormap, kw...)
    m = max(maximum(observed), maximum(modelled))
    lines!(ax, [0, m], [0, m]; color = :red, linewidth = 1)
    return ax, sc
end

# Short text in the top left corner of an axis, for the statistics of a panel
function _stats_note!(ax, lines...)
    text!(ax, 0.04, 0.96; text = join(lines, "\n"), space = :relative,
        align = (:left, :top), fontsize = 8)
end

"""
    plot_fit_to_observations(results, glacier; B = results.B, title = nothing) -> Figure

How well the model fits the observations, in one figure:

  - the ice thickness at each survey: model against observation, and the misfit divided by the
    uncertainty of each cell when the surveys have one;
  - the surface speed of the velocity map: model against observation (coloured by elevation),
    the mean misfit for each elevation band, and the map of the misfit;
  - the surface at each dated DEM: the map of model minus DEM, and the histograms.

Each part is drawn if the data exist. `B` is the bed of the model, as in
[`plot_surface_change`](@ref). Cells without data are light grey in the maps. For the model,
the observation and the difference at each date of one variable, see
[`plot_glacier_vs_observations`](@ref).
"""
function plot_fit_to_observations(
        results::Results, glacier;
        B::AbstractMatrix = results.B,
        title::Union{Nothing, AbstractString} = nothing
)
    td, sd = glacier.thicknessData, glacier.surfaceData
    has_H = !isnothing(td) && !isnothing(td.t) && !isempty(td.t)
    vd = glacier.velocityData
    has_V = !isempty(results.V) && !isnothing(vd) && !isempty(vd.vabs)
    outline = results.H[begin] .> 0
    nearest(t) = argmin(abs.(results.t .- t))
    return with_theme(paper_theme()) do
        fig = Figure(size = (1000, 1010))
        isnothing(title) || Label(fig[0, 1:3], title; fontsize = 12, font = :bold)
        row = 0
        letter = Iterators.Stateful('a':'z')
        lettered(s) = "($(popfirst!(letter))) $(s)"
        # Thickness
        if has_H
            row += 1
            all_chi = Float64[]
            all_err = Float64[]
            for (k, (t, Hobs)) in enumerate(zip(td.t, td.H))
                k > 2 && break
                H = results.H[nearest(t)]
                cells = Hobs .> 0
                r, b = _rmse_bias(H[cells] .- Hobs[cells])
                ax,
                _ = _scatter_fit!(fig[row, k], Hobs[cells], H[cells],
                    lettered("Thickness survey, $(_date_label(t))");
                    xlabel = "observed thickness (m)", ylabel = "modelled thickness (m)")
                _stats_note!(ax, "$(count(cells)) cells", "RMSE $(round(r; digits = 1)) m",
                    "bias $(round(b; digits = 1)) m")
                σ = isnothing(td.σ) ? fill(NaN, size(Hobs)) : td.σ[k]
                append!(all_err, H[cells] .- Hobs[cells])
                append!(all_chi, ((H .- Hobs) ./ σ)[cells .& isfinite.(σ)])
            end
            ax = Axis(fig[row, 3]; ylabel = "cells")
            if !isempty(all_chi)
                ax.title = lettered("Thickness misfit / uncertainty")
                ax.xlabel = "(model − observation) / σ"
                hist!(ax, all_chi; bins = 40, color = (:steelblue, 0.75))
                _stats_note!(ax, "mean χ² $(round(mean(abs2, all_chi); digits = 1))")
            else
                ax.title = lettered("Thickness misfit")
                ax.xlabel = "model − observation (m)"
                hist!(ax, all_err; bins = 40, color = (:steelblue, 0.75))
            end
            vlines!(ax, 0; color = :red, linewidth = 1)
        end
        # Velocity
        if has_V
            row += 1
            date1, date2 = first(vd.date1), first(vd.date2)
            t1, t2 = datetime_to_floatyear(date1), datetime_to_floatyear(date2)
            w = findall(t -> t1 <= t <= t2, results.t)
            isempty(w) && (w = [nearest((t1 + t2) / 2)])
            V̄ = sum(results.V[w]) ./ length(w)
            Vref = first(vd.vabs)
            cells = (Vref .> 0) .& outline
            r, b = _rmse_bias(V̄[cells] .- Vref[cells])
            gl = GridLayout()
            fig[row, 1] = gl
            ax,
            sc = _scatter_fit!(gl[1, 1], Vref[cells], V̄[cells],
                lettered("Surface speed, $(_date_label(date1)) – $(_date_label(date2))");
                xlabel = "observed speed (m/yr)", ylabel = "modelled speed (m/yr)",
                color = glacier.S[cells], colormap = :viridis)
            _stats_note!(ax, "$(count(cells)) cells", "RMSE $(round(r; digits = 1)) m/yr",
                "bias $(round(b; digits = 1)) m/yr")
            Colorbar(gl[1, 2], sc; label = "surface elevation (m)")
            # Mean misfit by elevation band
            z = glacier.S[cells]
            d = V̄[cells] .- Vref[cells]
            edges = range(floor(minimum(z) / 250) * 250, ceil(maximum(z) / 250) * 250; step = 250)
            mid, mean_d, std_d = Float64[], Float64[], Float64[]
            for i in 1:(length(edges) - 1)
                inband = (z .>= edges[i]) .& (z .< edges[i + 1])
                count(inband) < 5 && continue
                push!(mid, (edges[i] + edges[i + 1]) / 2)
                push!(mean_d, mean(d[inband]))
                push!(std_d, std(d[inband]))
            end
            ax2 = Axis(fig[row, 2]; title = lettered("Speed misfit by elevation band"),
                xlabel = "surface elevation (m)", ylabel = "model − observation (m/yr)")
            hlines!(ax2, 0; color = :red, linewidth = 1)
            errorbars!(ax2, mid, mean_d, std_d; whiskerwidth = 4, color = :gray40)
            scatter!(ax2, mid, mean_d; color = :black, markersize = 5)
            diff_map = ifelse.(cells, V̄ .- Vref, NaN)
            _map_panel!(fig[row, 3],
                MapPanel(diff_map, lettered("Speed: model − observation");
                    cmap = Reverse(:balance), label = "m/yr", symmetric = true,
                    nodata = BitMatrix(outline .& .!cells)),
                results,
                outline)
        end
        # Surface
        if !isnothing(sd)
            dems = findall(t -> results.t[begin] <= t <= results.t[end], sd.t)
            if !isempty(dems)
                row += 1
                free = _free_cells(glacier, B)
                labels = String[]
                residuals = Vector{Vector{Float64}}()
                for (c, k) in enumerate(dems)
                    c > 2 && break
                    ok = isfinite.(sd.S[k]) .& free
                    res_k = ifelse.(ok, _model_surface(results, B, sd.t[k]) .- sd.S[k], NaN)
                    r, b = _rmse_bias(res_k[ok])
                    _map_panel!(fig[row, c],
                        MapPanel(res_k,
                            lettered("Surface, model − DEM, $(_date_label(sd.t[k]))\n" *
                                     "RMSE $(round(r; digits = 1)) m, bias $(round(b; digits = 1)) m");
                            cmap = Reverse(:balance), label = "m", symmetric = true,
                            nodata = BitMatrix(free .& .!ok)),
                        results,
                        outline)
                    push!(residuals, res_k[ok])
                    push!(labels, _date_label(sd.t[k]))
                end
                ax3 = Axis(fig[row, 3]; title = lettered("Surface residuals"),
                    xlabel = "model − DEM (m)", ylabel = "cells")
                for (res_k, lab, col) in zip(residuals, labels, (:steelblue, :darkorange))
                    hist!(ax3, res_k; bins = 50, color = (col, 0.6), label = lab)
                end
                vlines!(ax3, 0; color = :red, linewidth = 1)
                axislegend(ax3; position = :rt)
            end
        end
        colgap!(fig.layout, 12)
        rowgap!(fig.layout, 12)
        fig
    end
end
