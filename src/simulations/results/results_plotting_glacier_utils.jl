###############################################################################
# Glacier-level plotting utilities (CairoMakie)
#
# Layout:
#   1. Orientation helper:   reverseForHeatmap
#   2. Internal helpers:     _resolve_temporal_snapshot, _results_plot_metadata,
#                            _overlay_contour!, _decorate_geo_axis!
#   3. Public plot functions: plot_glacier_heatmaps, plot_glacier_quivers,
#                             plot_glacier_difference_evolution,
#                             plot_glacier_statistics_evolution,
#                             plot_glacier_integrated_volume, plot_bias
#   4. Router:               plot_glacier (dispatches on plot_type string/Symbol)
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 1. Orientation helper
# ─────────────────────────────────────────────────────────────────────────────

"""
    reverseForHeatmap(inp, x, y) -> Matrix

Out-of-place reverse of a matrix so that the heatmap is displayed in the correct
geographic orientation. Flips along each axis whose coordinate vector is descending.

# Arguments

  - `inp::Matrix{F}`: The matrix to reverse.
  - `x::Vector{F}`: Values of the x axis.
  - `y::Vector{F}`: Values of the y axis.

# Returns

  - Out-of-place copy of inp that has been reversed if needed.
"""
function reverseForHeatmap(
        inp::Matrix{F},
        x::Vector{F},
        y::Vector{F}
) where {F <: AbstractFloat}
    reverseX = x[end]<x[begin]
    reverseY = y[end]<y[begin]
    data = copy(inp)
    if reverseX
        data = reverse(data, dims = 1)
    end
    if reverseY
        data = reverse(data, dims = 2)
    end
    return data
end

# ─────────────────────────────────────────────────────────────────────────────
# 2. Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    _resolve_temporal_snapshot(data, timeIdx) -> Matrix

If `data` is a vector of matrices, pick element `timeIdx` (or the last one when
`timeIdx` is `nothing`). Otherwise return `data` as-is.
"""
function _resolve_temporal_snapshot(
        data,
        timeIdx::Union{Nothing, Int64};
        varname::String = "Data"
)
    if data isa AbstractVector
        @assert !isempty(data) "$(varname) is an empty vector"
        @assert isnothing(timeIdx) || length(data) >= timeIdx "The provided index=$(timeIdx) is greater than the size of the vector for $(varname) which is $(length(data))"
        return isnothing(timeIdx) ? data[end] : data[timeIdx]
    end
    return data
end

"""
    _results_plot_metadata(results; plotContour) -> NamedTuple

Extract common metadata (lon, lat, x, y, rgi_id, Δx, mask, contour) from a
`Results` object. The glacier mask is derived from the initial thickness field.
"""
function _results_plot_metadata(results::Results; plotContour::Bool = false)
    mask = results.H[begin] .> 0.0
    nx, ny = size(results.H[begin])
    ctr = plotContour ? Contour.contour(collect(1:nx), 1+ny .- collect(1:ny), mask, 0.5) :
          nothing
    return (
        lon = results.lon,
        lat = results.lat,
        x = results.x,
        y = results.y,
        rgi_id = results.rgi_id,
        Δx = results.Δx,
        mask = mask,
        contour = ctr
    )
end

"""
    _overlay_contour!(ax, ctr)

Draw glacier-boundary contour lines on an existing axis. No-op when `ctr` is `nothing`.
"""
function _overlay_contour!(ax::Axis, ctr)
    isnothing(ctr) && return
    for curve in ctr.lines
        xs = first.(curve.vertices)
        ys = last.(curve.vertices)
        lines!(ax, xs, ys, color = :black, linewidth = 1)
    end
end

"""
    _decorate_geo_axis!(ax, nx, ny, lon, lat, Δx; scale_text_size, num_vars)

Add geographic labels (central lon/lat tick), axis styling, and a Δx-based
scale bar to a glacier heatmap axis. Shared by heatmap, difference-evolution,
and gridded-data plotting functions.
"""
function _decorate_geo_axis!(ax::Axis, nx, ny, lon, lat, Δx;
        scale_text_size::Union{Nothing, Float64} = nothing,
        num_vars::Int = 1)
    ax.xlabel = "Longitude"
    ax.ylabel = "Latitude"
    ax.xticks = ([round(nx / 2)], ["$(round(lon; digits=6)) °"])
    ax.yticks = ([round(ny / 2)], ["$(round(lat; digits=6)) °"])
    ax.yticklabelrotation = π / 2
    ax.ylabelpadding = 5
    ax.yticklabelalign = (:center, :bottom)

    # Scale bar
    scale_width = 0.10 * nx
    scale_number = round(Δx * scale_width / 1000; digits = 1)
    if isnothing(scale_text_size)
        textsize = if num_vars <= 1
            1.2 * scale_width
        elseif num_vars == 2
            0.9 * scale_width
        else
            0.5 * scale_width
        end
    else
        textsize = scale_text_size
    end
    poly!(
        ax, Rect(nx - round(0.15 * nx), round(0.075 * ny), scale_width, scale_width / 10),
        color = :black)
    text!(ax, "$scale_number km",
        position = (nx - round(0.15 * nx) + scale_width / 16,
            round(0.075 * ny) + scale_width / 10),
        fontsize = textsize)
end

function _fixed_colorbar_ticks(lower::Real, upper::Real; count::Int = 5)
    if !isfinite(lower) || !isfinite(upper)
        return [0.0, 1.0]
    end
    if upper <= lower
        return [float(lower)]
    end
    return collect(range(float(lower), float(upper), length = count))
end

# ─────────────────────────────────────────────────────────────────────────────
# 3. Public plot functions
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_glacier_heatmaps(
        results::Results,
        variables::Vector{Symbol},
        title_mapping::Dict;
        scale_text_size::Union{Nothing,Float64}=nothing,
        timeIdx::Union{Nothing,Int64}=nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        plotContour::Bool=false,
    ) -> Figure

Plot heatmaps for glacier variables.

# Arguments

  - `results::Results`: The results object containing the data to be plotted.
  - `variables::Vector{Symbol}`: A list of variables to be plotted.
  - `title_mapping::Dict`: A dictionary mapping variable names to their titles and colormaps.
  - `scale_text_size::Union{Nothing,Float64}`: Optional argument to scale the text size.
  - `timeIdx::Union{Nothing,Int64}`:: Optional argument to select the index at which
    data should be plotted when dealing with vector of matrix. Default is nothing
    which selects the last element available.
  - `figsize::Union{Nothing, Tuple{Int64, Int64}}`: Size of the figure.
  - `plotContour::Bool`: Whether to add a contour plot representing the glacier borders at
    the beginning of the simulation on top of each of the figures. Default is false.

# Returns

  - A plot of the glacier heatmaps.
"""
function plot_glacier_heatmaps(
        results::Results,
        variables::Vector{Symbol},
        title_mapping::Dict;
        scale_text_size::Union{Nothing, Float64} = nothing,
        timeIdx::Union{Nothing, Int64} = nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        plotContour::Bool = false
)
    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)

    # Dictionary of variable-specific colormaps
    colormap_mapping = Dict(key => value[3] for (key, value) in title_mapping)

    # Extract metadata about the glacier
    meta = _results_plot_metadata(results; plotContour = plotContour)
    lon, lat = meta.lon, meta.lat
    x, y = meta.x, meta.y
    rgi_id, Δx = meta.rgi_id, meta.Δx
    mask = meta.mask

    ice_thickness_vars = [:H, :H₀, :H_glathida, :H_ref] # Ice thickness variables
    velocity_vars = [:V, :Vx, :Vy, :V_ref] # Velocity variables, excluding V_diff

    # Initialize max_values for ice thickness and velocity separately, considering only given variables
    max_values_ice = []
    max_values_velocity = []

    for var in intersect(union(ice_thickness_vars, velocity_vars), variables)  # Check only given vars for maximum
        if hasproperty(results, var)
            current_matrix = getfield(results, var)
            if !isnothing(current_matrix) && !isempty(current_matrix)
                current_snapshot = _resolve_temporal_snapshot(current_matrix, timeIdx; varname = string(var))
                current_snapshot = copy(current_snapshot)
                current_snapshot[.!mask] .= NaN
                finite_vals = filter(isfinite, vec(current_snapshot))
                maxval = isempty(finite_vals) ? 0.0 : maximum(finite_vals)
                if var in ice_thickness_vars
                    push!(max_values_ice, maxval)
                elseif var in velocity_vars
                    push!(max_values_velocity, maxval)
                end
            end
        end
    end

    # Determine global maximum for ice and velocity separately.
    # Guard against zero or NaN maxima (e.g. fully masked snapshots) so that
    # CairoMakie never receives a degenerate colorrange like (0, 0).
    global_max_ice = isempty(max_values_ice) ? nothing :
                     max(maximum(max_values_ice), eps(Float64))
    global_max_velocity = isempty(max_values_velocity) ? nothing :
                          max(maximum(max_values_velocity), eps(Float64))

    num_vars = length(variables)
    rows, cols = if num_vars == 1
        2, 2
    elseif num_vars == 2
        3, 2
    elseif num_vars in [3, 4]
        3, 4
    else
        error("Unsupported number of variables.")
    end

    figKwargs[:layout] = GridLayout(rows, cols)
    fig = Figure(; figKwargs...)
    for (i, var) in enumerate(variables)
        ax_row = div(i - 1, 2) + 1
        ax_col = 2 * (rem(i - 1, 2)) + 1
        ax = Axis(fig[ax_row, ax_col], aspect = DataAspect())
        data = _resolve_temporal_snapshot(getfield(results, var), timeIdx; varname = string(var))
        title, unit = get(title_mapping, string(var), (string(var), ""))
        data = copy(data)

        nx, ny = size(data)
        colormap = get(colormap_mapping, string(var), :cool)  # Default colormap

        if (var in ice_thickness_vars) || (var in velocity_vars)
            data[.!mask] .= NaN
        end
        if var==:H_glathida
            # For GlaThiDa variable, replace zeros by NaN
            data[mask .& (data .== 0)] .= NaN
        end

        # Apply global_max_ice to ice thickness variables and global_max_velocity to velocity variables
        cb_ticks = nothing
        if var in ice_thickness_vars
            colorrange = (0.0, global_max_ice)
            hm = heatmap!(ax, reverseForHeatmap(data, x, y),
                colormap = colormap, colorrange = colorrange)
            cb_ticks = _fixed_colorbar_ticks(colorrange[1], colorrange[2])
        elseif var in velocity_vars
            colorrange = (0.0, global_max_velocity)
            hm = heatmap!(ax, reverseForHeatmap(data, x, y), colormap = colormap,
                colorrange = colorrange)
            cb_ticks = _fixed_colorbar_ticks(colorrange[1], colorrange[2])
        else
            hm = heatmap!(ax, reverseForHeatmap(data, x, y), colormap = colormap)
            finite_vals = filter(isfinite, vec(data))
            if !isempty(finite_vals)
                lower, upper = extrema(finite_vals)
                cb_ticks = _fixed_colorbar_ticks(lower, upper)
            end
        end
        cb_kwargs = isnothing(cb_ticks) ? NamedTuple() : (; ticks = cb_ticks)
        cb = Colorbar(fig[ax_row, ax_col + 1], hm; cb_kwargs...)
        Observables.connect!(
            cb.height, @lift CairoMakie.Fixed($(viewport(ax.scene)).widths[2]))
        Label(fig[ax_row, ax_col + 1], "$var ($unit)",
            fontsize = 14, valign = :top, padding = (0, -25))

        _overlay_contour!(ax, meta.contour)
        ax.title = "$title"
        _decorate_geo_axis!(ax, nx, ny, lon, lat, Δx;
            scale_text_size = scale_text_size, num_vars = num_vars)
    end

    fig[0, :] = Label(fig, "$rgi_id")
    resize_to_layout!(fig)
    return fig
end

"""
    plot_glacier_quivers(
        results::Results,
        variables::Vector{Symbol},
        title_mapping::Dict;
        timeIdx::Union{Nothing,Int64} = nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        lengthscale::Float64 = 0.00001,
        tiplength::Float64 = 0.5,
    ) -> Figure

Plot quivers for glacier variables.

# Arguments

  - `results::Results`: The results object containing the data to be plotted.
  - `variables::Vector{Symbol}`: A list of variables to be plotted.
  - `title_mapping::Dict`: A dictionary mapping variable names to their titles and colormaps.
  - `timeIdx::Union{Nothing,Int64}`:: Optional argument to select the index at which
    data should be plotted when dealing with vector of matrix. Default is nothing
    which selects the last element available.
  - `figsize::Union{Nothing, Tuple{Int64, Int64}}`: Size of the figure.
  - `lengthscale::Float64`: Lengthscale of the arrows in the quiver plot.
  - `tiplength::Float64`: Length of the arrow in the quiver plot.

# Returns

  - A plot of the glacier quivers.
"""
function plot_glacier_quivers(
        results::Results,
        variables::Vector{Symbol},
        title_mapping::Dict;
        timeIdx::Union{Nothing, Int64} = nothing,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
        lengthscale::Float64 = 0.00001,
        tiplength::Float64 = 0.5
)
    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)

    # Extract metadata about the glacier
    meta = _results_plot_metadata(results; plotContour = false)
    x, y = meta.x, meta.y

    num_vars = length(variables)
    rows, cols = if num_vars == 1
        1, 1
    elseif num_vars == 2
        1, 2
    else
        error("Unsupported number of variables.")
    end

    figKwargs[:layout] = GridLayout(rows, cols)
    fig = Figure(; figKwargs...)
    for (ax_col, var) in enumerate(variables)
        ax = Axis(fig[1, ax_col], aspect = DataAspect())
        title = get(title_mapping, string(var), (string(var), ""))[1]

        Vx = _resolve_temporal_snapshot(getfield(results, var == :V ? :Vx : :Vx_ref),
            timeIdx; varname = string(var) * "x")
        Vy = _resolve_temporal_snapshot(getfield(results, var == :V ? :Vy : :Vy_ref),
            timeIdx; varname = string(var) * "y")

        X, Y = meshgrid(x, y)

        positions = Point2f.(reshape(X, :), reshape(Y, :))
        directions = Vec2f.(Vx, -Vy)
        arrows2d!(
            ax, positions, directions; tiplength = tiplength, lengthscale = lengthscale)

        ax.title = "$title"
        ax.xlabel = "Longitude"
        ax.ylabel = "Latitude"
        ax.yticklabelrotation = π/2
        ax.ylabelpadding = 5
        ax.yticklabelalign = (:center, :bottom)
    end

    resize_to_layout!(fig)
    return fig
end

"""
    plot_glacier_difference_evolution(
        results::Results,
        variables::Vector{Symbol},
        title_mapping;
        tspan::Tuple{F,F}=results.tspan,
        metrics::Vector{String}="difference",
        figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing,
    ) where {F<:AbstractFloat}

Plot the evolution of the difference in a glacier variable over time.

# Arguments

  - `results::Results`: The simulation results object containing the data to be plotted.
  - `variables::Vector{Symbol}`: The variable to be plotted.
  - `title_mapping`: A dictionary mapping variable names to their titles.
  - `tspan::Tuple{F,F}`: A tuple representing the start and end time for the simulation.
  - `metrics::Vector{String}`: Metrics to visualize, e.g., `["difference"]`.
  - `figsize::Union{Nothing, Tuple{Int64, Int64}}`: Size of the figure.

# Returns

  - A plot of the glacier difference evolution.
"""
function plot_glacier_difference_evolution(
        results::Results,
        variables::Vector{Symbol},
        title_mapping;
        tspan::Tuple{F, F} = results.tspan,
        metrics::Vector{String} = "difference",
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing
) where {F <: AbstractFloat}
    # Check if more than one variable is passed
    @assert length(variables) == 1 "Only one variable can be passed to this function."

    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)

    # Check for valid metrics
    valid_metrics = ["hist", "difference"]
    for metric in metrics
        if !(metric in valid_metrics)
            error("Invalid metric: $metric. Valid metrics are: $valid_metrics")
            return
        end
    end

    # Extract data for the variable
    data = getfield(results, variables[1])

    # Extract metadata about the glacier
    lon = results.lon
    lat = results.lat
    x = results.x
    y = results.y
    rgi_id = results.rgi_id
    Δx = results.Δx

    # Check the shape of the extracted data
    @assert data isa AbstractVector{<:AbstractMatrix{<:AbstractFloat}} "Only temporal quantities can be used in this function."

    # Print plot information
    variable_title = get(title_mapping, variables[1], variables[1])
    data_diff = data[end] - data[begin]
    mask = results.H[begin] .> 0.0
    data_diff[.!mask] .= NaN

    # Determine whether to create a single plot or a subplot
    if metrics == ["hist"]
        fig = Figure(; figKwargs...)
        ax = Axis(
            fig[1, 1],
            xlabel = "Δ$variable_title ($(title_mapping[string(variables[1])][2]))",
            ylabel = "Frequency",
            title = "Histogram of $(title_mapping[string(variables[1])][1]) evolution ($rgi_id)"
        )
    elseif metrics == ["difference"]
        fig = Figure(; figKwargs...)
        ax_diff = Axis(fig[1, 1],
            title = "$(title_mapping[string(variables[1])][1]) evolution ($rgi_id)",
            aspect = DataAspect())
    else
        figKwargs[:layout] = GridLayout(1, 3)
        fig = Figure(; figKwargs...)
        ax = Axis(
            fig[1, 3],
            xlabel = "Δ$variable_title ($(title_mapping[string(variables[1])][2]))",
            ylabel = "Frequency",
            title = "Histogram of $(title_mapping[string(variables[1])][1])\n evolution"
        )
        ax_diff = Axis(
            fig[1, 1], title = "$(title_mapping[string(variables[1])][1]) evolution",
            aspect = DataAspect())
        fig[0, :] = Label(fig, "$rgi_id")
    end

    # Plot based on the metric
    for metric in metrics
        if metric == "hist"
            hist!(ax, vec(data_diff[mask]), bins = 50)
            ax.limits[] = (
                minimum(data_diff[mask]), maximum(data_diff[mask]), nothing, nothing)
        elseif metric == "difference"
            nx, ny = size(data_diff)

            # Calculate the symmetric color range
            max_abs_value = max(
                abs(minimum(data_diff[mask])), abs(maximum(data_diff[mask])))

            hm_diff = heatmap!(
                ax_diff, reverseForHeatmap(data_diff, x, y), colormap = :redsblues,
                colorrange = (-max_abs_value, max_abs_value))

            _decorate_geo_axis!(ax_diff, nx, ny, lon, lat, Δx;
                num_vars = length(metrics))

            cb = Colorbar(fig[1, 2], hm_diff)
            Observables.connect!(
                cb.height, @lift CairoMakie.Fixed($(viewport(ax_diff.scene)).widths[2]))
            Label(
                fig[1, 2], "Δ$variable_title ($(title_mapping[string(variables[1])][2]))",
                fontsize = 14, valign = :top, padding = (0, -25))
        end
    end

    resize_to_layout!(fig)
    return fig  # Return the main figure
end

"""
    plot_glacier_statistics_evolution(
        results::Results,
        variables::Vector{Symbol},
        title_mapping;
        metrics="median",
        tspan,
        threshold=0.5,
        figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing,
    )

Plot the evolution of statistics for multiple glacier variables over time.

# Arguments

  - `results::Results`: The simulation results object containing the data to be plotted.
  - `variables::Vector{Symbol}`: A list of variables to be plotted.
  - `title_mapping`: A dictionary mapping variable names to their titles.
  - `metrics`: Metrics to visualize, e.g., "average", "median", "min", "max", and "std". Default is "median".
  - `tspan`: A tuple representing the start and end time for the simulation.
  - `threshold`: A threshold value to filter the data. Default is 0.5.
  - `figsize::Union{Nothing, Tuple{Int64, Int64}}`: Size of the figure.

# Returns

  - A plot of the glacier statistics evolution.
"""
function plot_glacier_statistics_evolution(
        results::Results,
        variables::Vector{Symbol},
        title_mapping;
        metrics = "median",
        tspan,
        threshold = 0.5,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing
)
    # Check if more than one variable is passed
    @assert length(variables) == 1 "Only one variable can be passed to this function."

    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)

    # Extract data for the variable
    data = getfield(results, variables[1])

    # Check the shape of the extracted data
    @assert data isa AbstractVector{<:AbstractMatrix{<:AbstractFloat}} "Only temporal quantities can be used in this function."

    # Check for valid metrics
    valid_metrics = ["average", "median", "max", "std", "min"]
    for metric in metrics
        if !(metric in valid_metrics)
            error("Invalid metric: $metric. Valid metrics are: $valid_metrics")
            return
        end
    end

    # Extract the rgi_id
    rgi_id = results.rgi_id

    # Create a time vector
    t = range(tspan[1], stop = tspan[2], length = length(getfield(results, variables[1])))

    # Create a single plot for all other metrics
    fig = Figure(; figKwargs...)
    ax = Axis(fig[1, 1],
        xlabel = "Time (years)",
        ylabel = "$(title_mapping[string(variables[1])][1]) ($(title_mapping[string(variables[1])][2]))",
        title = "Metrics for $(title_mapping[string(variables[1])][1]) through time ($rgi_id)")

    # If "average" or "std" is in metrics, calculate them
    if "average" in metrics || "std" in metrics
        avg_vals = [let v = filter(x -> !isnan(x) && x >= threshold, matrix[:])
                        isempty(v) ? NaN : mean(v)
                    end
                    for matrix in data]
        std_vals = [let v = filter(x -> !isnan(x) && x >= threshold, matrix[:])
                        isempty(v) ? NaN : std(v)
                    end
                    for matrix in data]
    end

    # Calculate and plot metrics
    for metric in metrics
        if metric == "average"
            if "std" in metrics
                band!(ax, t, avg_vals .- std_vals, avg_vals .+ std_vals,
                    alpha = 0.1, label = "Std Dev", color = :lightgray)
            end
            lines!(ax, t, avg_vals, label = "Average")
        elseif metric == "median"
            median_vals = [let v = filter(x -> !isnan(x) && x >= threshold, matrix[:])
                               isempty(v) ? NaN : median(v)
                           end
                           for matrix in data]
            lines!(ax, t, median_vals, label = "Median")
        elseif metric == "min"
            min_vals = [let v = filter(x -> !isnan(x) && x >= threshold, matrix[:])
                            isempty(v) ? NaN : minimum(v)
                        end
                        for matrix in data]
            lines!(ax, t, min_vals, linestyle = :dot, label = "Min")
        elseif metric == "max"
            max_vals = [let v = filter(x -> !isnan(x) && x >= threshold, matrix[:])
                            isempty(v) ? NaN : maximum(v)
                        end
                        for matrix in data]
            lines!(ax, t, max_vals, linestyle = :dot, label = "Max")
        end
    end

    leg = Legend(fig, ax)
    fig[1, 2] = leg
    resize_to_layout!(fig)

    return fig  # Return the main figure
end

"""
    plot_glacier_integrated_volume(
        results,
        variables,
        title_mapping;
        tspan,
        figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing,
    )

Plot the integrated volume of a glacier variable over time.

# Arguments

  - `results::Results`: The results object containing the data to be plotted.
  - `variables::Vector{Symbol}`: The variable to be plotted.
  - `title_mapping`: A dictionary mapping variable names to their titles.
  - `tspan`: A tuple representing the start and end time for the simulation.
  - `figsize::Union{Nothing, Tuple{Int64, Int64}}`: Size of the figure.

# Returns

  - A plot of the glacier integrated volume.
"""
function plot_glacier_integrated_volume(
        results,
        variables,
        title_mapping;
        tspan,
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing
)

    # Determine pixel area
    area=results.Δx*results.Δy

    # Check if more than one variable is passed
    @assert length(variables) == 1 "Only one variable can be passed to this function."

    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)

    # Extract the rgi_id
    rgi_id = results.rgi_id

    # Create a time vector
    t = range(tspan[1], stop = tspan[2], length = length(getfield(results, variables[1])))

    # Extract data for the variable
    data = getfield(results, variables[1])

    # Check the shape of the extracted data
    @assert data isa AbstractVector{<:AbstractMatrix{<:AbstractFloat}} "Only temporal quantities can be used in this function."

    # Calculate integrated ice volume for each time step
    integrated_ice_volume = [sum(matrix) * area for matrix in data]/(10^9)  # Multiply by area and convert to km^3

    nDigitsYticks = max(2,
        Int(round(-log10(maximum(integrated_ice_volume) - minimum(integrated_ice_volume)))))

    # Plot the integrated ice volume as a function of time
    fig = Figure(; figKwargs...)
    ax = Axis(fig[1, 1]; ytickformat = x -> string.(round.(x; digits = nDigitsYticks)))
    lines!(ax, t, integrated_ice_volume, color = :blue)
    ax.xlabel = "Time (years)"
    ax.ylabel = "Integrated ice volume (km³)   "
    ax.title = "Evolution of integrated ice volume ($rgi_id)"

    resize_to_layout!(fig)
    return fig  # Return the main figure with the plot
end

"""
    plot_bias(
        results,
        variables;
        treshold = [0, 0],
        figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing,
    )

Plot the bias of the glacier integrated volume over the specified time span.

# Arguments

  - `results::Results`: The results object containing the data to be plotted.
  - `variables::Vector{Symbol}`: The variables to be plotted.
  - `title_mapping::Dict{Symbol, String}`: A dictionary mapping variable names to their titles.
  - `tspan::Tuple{Float64, Float64}`: A tuple representing the start and end time for the simulation.
  - `figsize::Union{Nothing, Tuple{Int64, Int64}}`: Size of the figure.

# Returns

  - A plot of the glacier integrated volume bias.
"""
function plot_bias(
        results,
        variables;
        treshold = [0, 0],
        figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing
)
    # Check for exactly two variables
    @assert length(variables) == 2 "Exactly two variables are required for the scatter plot."

    figKwargs = isnothing(figsize) ? Dict{Symbol, Any}() :
                Dict{Symbol, Any}(:size => figsize)

    # Ensure treshold is an array of length 2
    if length(treshold) == 1
        treshold = [treshold[1], treshold[1]]
    end

    # Extract data from results
    rgi_id = results.rgi_id
    x_values = getfield(results, variables[1])
    y_values = getfield(results, variables[2])

    # Filter non-zero observations if necessary
    if :H_obs in variables || :V_obs in variables
        obs_key = :H_obs in variables ? :H_obs : :V_obs
        mask_H = getfield(results, obs_key) .>
                 treshold[1] * maximum(getfield(results, obs_key))
        mask_H_2 = getfield(results, obs_key) .<
                   (1 - treshold[2]) * maximum(getfield(results, obs_key))
        x_values = x_values[mask_H .& mask_H_2]
        y_values = y_values[mask_H .& mask_H_2]
    end

    # Calculate metrics
    differences = x_values .- y_values
    rmse = sqrt(mean(differences .^ 2)) # Root Mean Square Error
    bias = mean(differences) # Bias
    ss_res = sum(differences .^ 2) # Sum of squares of residuals
    ss_tot = sum((x_values .- mean(x_values)) .^ 2) # Total sum of squares
    r_squared = 1 - (ss_res / ss_tot) # R-squared

    # Plotting
    fig = Figure(; figKwargs...)
    ax = Axis(fig[1, 1], xlabel = string(variables[1]), ylabel = string(variables[2]),
        title = "Scatter plot for RGI ID: " * rgi_id)

    scatter!(
        ax, vec(x_values), vec(y_values), markersize = 5, color = :blue, label = "Data")
    xmin, xmax = minimum(x_values), maximum(x_values)
    ymin, ymax = minimum(y_values), maximum(y_values)
    lines!(ax, [xmin, xmax], [xmin, xmax], linestyle = :dash, color = :red, label = "y = x")

    # Display metrics on the plot
    metrics_text = "RMSE: $(round(rmse, digits=2))\nR²: $(round(r_squared, digits=2))\nBias: $(round(bias, digits=2))"
    text!(ax, metrics_text, position = (xmax, ymax), align = (:right, :top), color = :black)

    return fig
end

# ─────────────────────────────────────────────────────────────────────────────
# 4. Router
# ─────────────────────────────────────────────────────────────────────────────

"""
    plot_glacier(results, plot_type, variables; kwargs...) -> Figure

High-level entry point that dispatches to specific glacier plotting functions
based on `plot_type`:

  - `"heatmaps"` → `plot_glacier_heatmaps`
  - `"quivers"` → `plot_glacier_quivers`
  - `"evolution difference"` → `plot_glacier_difference_evolution`
  - `"evolution statistics"` → `plot_glacier_statistics_evolution`
  - `"integrated volume"` → `plot_glacier_integrated_volume`
  - `"bias"` → `plot_bias`
  - `"dem"` → `plot_glacier_dem`
"""
function plot_glacier(
        results::Results,
        plot_type::Union{String, Symbol},
        variables::Vector{Symbol};
        kwargs...)
    title_mapping = Dict(
        "H" => ("Predicted ice thickness", "m", :YlGnBu),
        "H₀" => ("Ice thickness", "m", :YlGnBu),
        "H_glathida" => ("Ice thickness (GlaThiDa)", "m", :YlGnBu),
        "S" => ("Surface topography", "m", :terrain),
        "B" => ("Bed topography", "m", :terrain),
        "Vx" => ("Ice surface velocity\n(X-direction)", "m/y", :viridis),
        "Vy" => ("Ice surface velocity\n(Y-direction)", "m/y", :viridis),
        "H_ref" => ("Observed ice thickness", "m", :YlGnBu),
        "H_diff" => ("Ice thickness difference", "m", :RdBu),
        "V" => ("Predicted ice surface velocity", "m/y", :viridis),
        "V_ref" => ("Observed ice surface velocity", "m/y", :viridis),
        "V_diff" => ("Ice surface velocity difference", "m/y", :RdBu)
    )

    plot_type_s = plot_type isa Symbol ? String(plot_type) : plot_type

    if plot_type_s == "heatmaps"
        return plot_glacier_heatmaps(results, variables, title_mapping; kwargs...)
    elseif plot_type_s == "quivers"
        return plot_glacier_quivers(results, variables, title_mapping; kwargs...)
    elseif plot_type_s == "evolution difference"
        return plot_glacier_difference_evolution(
            results, variables, title_mapping; kwargs...)
    elseif plot_type_s == "evolution statistics"
        return plot_glacier_statistics_evolution(
            results, variables, title_mapping; kwargs...)
    elseif plot_type_s == "integrated volume"
        return plot_glacier_integrated_volume(results, variables, title_mapping; kwargs...)
    elseif plot_type_s == "bias"
        return plot_bias(results, variables; kwargs...)
    elseif plot_type_s == "dem"
        return plot_glacier_dem(results; kwargs...)
    else
        error("Invalid plot_type: $plot_type_s")
    end
end
