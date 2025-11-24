export plot_glacier, plot_glacier_heatmaps, plot_glacier_quivers, plot_glacier_difference_evolution, plot_glacier_statistics_evolution, plot_glacier_integrated_volume
export plot_gridded_data

using CairoMakie: Axis


"""
    reverseForHeatmap(
        inp::Matrix{F},
        x::Vector{F},
        y::Vector{F}
    ) where {F <: AbstractFloat}

Out-of-place reverse of a matrix based on the values of the x and y axes.
This function corrects the orientation so that the heatmap is displayed correctly.

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
        data = reverse(data, dims=1)
    end
    if reverseY
        data = reverse(data, dims=2)
    end
    return data
end

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
    scale_text_size::Union{Nothing,Float64}=nothing,
    timeIdx::Union{Nothing,Int64}=nothing,
    figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
    plotContour::Bool=false,
)
    figKwargs = isnothing(figsize) ? Dict{Symbol,Any}() : Dict{Symbol,Any}(:size => figsize)

    # Dictionary of variable-specific colormaps
    colormap_mapping = Dict(key => value[3] for (key, value) in title_mapping)

    # Extract metadata about the glacier
    lon = results.lon
    lat = results.lat
    x = results.x
    y = results.y
    rgi_id = results.rgi_id
    Δx = results.Δx
    mask = results.H[begin] .> 0.0
    nx, ny = size(results.H[begin])
    ctr = plotContour ? Contour.contour(collect(1:nx),1+ny.-collect(1:ny),mask, 0.5) : nothing

    ice_thickness_vars = [:H, :H₀, :H_glathida, :H_ref] # Ice thickness variables
    velocity_vars = [:V, :Vx, :Vy, :V_ref] # Velocity variables, excluding V_diff

    # Initialize max_values for ice thickness and velocity separately, considering only given variables
    max_values_ice = []
    max_values_velocity = []

    for var in intersect(union(ice_thickness_vars, velocity_vars), variables)  # Check only given vars for maximum
        if hasproperty(results, var)
            current_matrix = getfield(results, var)
            if !isnothing(current_matrix) && !isempty(current_matrix)
                if typeof(current_matrix) <: Vector
                    @assert length(current_matrix)>0 "Variable $(var) is an empty vector"
                    @assert (isnothing(timeIdx)) || (size(current_matrix,1)>=timeIdx) "The provided index=$(timeIdx) is greater than the size of the vector for $(var) which is $(size(current_matrix,1))"
                    maxval = maximum(replace(isnothing(timeIdx) ? current_matrix[end] : current_matrix[timeIdx], NaN => 0.))
                    # maxval = maximum([maximum(replace(M, NaN => 0.)) for M in current_matrix])
                else
                    maxval = maximum(replace(current_matrix, NaN => 0.))
                end
                if var in ice_thickness_vars
                    push!(max_values_ice, maxval)
                elseif var in velocity_vars
                    push!(max_values_velocity, maxval)
                end
            end
        end
    end

    # Determine global maximum for ice and velocity separately
    global_max_ice = isempty(max_values_ice) ? nothing : maximum(max_values_ice)
    global_max_velocity = isempty(max_values_velocity) ? nothing : maximum(max_values_velocity)

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
        ax = Axis(fig[ax_row, ax_col], aspect=DataAspect())
        data = deepcopy(getfield(results, var))
        title, unit = get(title_mapping, string(var), (string(var), ""))

        if typeof(data) <: Vector
            @assert length(data)>0 "Variable $(var) is an empty vector"
            @assert (isnothing(timeIdx)) || (size(data,1)>=timeIdx) "The provided index=$(timeIdx) is greater than the size of the vector for $(var) which is $(size(data,1))"
            data = isnothing(timeIdx) ? data[end] : data[timeIdx]
        end

        nx, ny = size(data)
        colormap = get(colormap_mapping, string(var), :cool)  # Default colormap

        mask = results.H[begin] .> 0.0
        if (var in ice_thickness_vars) || (var in velocity_vars)
            data[.!mask] .= NaN
        end
        if var==:H_glathida
            # For GlaThiDa variable, replace zeros by NaN
            data[mask .& (data.==0)] .= NaN
        end

        # Apply global_max_ice to ice thickness variables and global_max_velocity to velocity variables
        if var in ice_thickness_vars
            hm = heatmap!(ax, reverseForHeatmap(data, x, y), colormap=colormap, colorrange=(0, global_max_ice))
        elseif var in velocity_vars
            hm = heatmap!(ax, reverseForHeatmap(data, x, y), colormap=colormap, colorrange=(0, global_max_velocity))
        else
            hm = heatmap!(ax, reverseForHeatmap(data, x, y), colormap=colormap)
        end
        cb = Colorbar(fig[ax_row, ax_col + 1], hm)
        Observables.connect!(cb.height, @lift CairoMakie.Fixed($(viewport(ax.scene)).widths[2]))
        Label(fig[ax_row, ax_col + 1], "$var ($unit)", fontsize=14, valign=:top, padding=(0, -25))

        if plotContour
            for curve in ctr.lines
                xs = first.(curve.vertices)
                ys = last.(curve.vertices)
                lines!(ax, xs, ys, color=:black, linewidth=1)
            end
        end

        ax.title = "$title"
        ax.xlabel = "Longitude"
        ax.ylabel = "Latitude"
        ax.xticks=([round(nx/2)], ["$(round(lon;digits=6)) °"])
        ax.yticks=([round(ny/2)], ["$(round(lat;digits=6)) °"])
        ax.yticklabelrotation = π/2
        ax.ylabelpadding = 5
        ax.yticklabelalign = (:center, :bottom)

        scale_width = 0.10*nx
        scale_number = round(Δx * scale_width / 1000; digits=1) # Convert to km
        if scale_text_size === nothing
            textsize = if num_vars == 1
                1.2*scale_width
            elseif num_vars == 2
                0.9*scale_width
            else
                0.5*scale_width
            end
        else
            textsize = scale_text_size
        end

        poly!(ax, Rect(nx - round(0.15*nx), round(0.075*ny), scale_width, scale_width/10), color=:black)
        text!(ax, "$scale_number km", position=(nx - round(0.15*nx) + scale_width/16, round(0.075*ny) + scale_width/10), fontsize=textsize)
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
    timeIdx::Union{Nothing,Int64} = nothing,
    figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
    lengthscale::Float64 = 0.00001,
    tiplength::Float64 = 0.5,
)
    figKwargs = isnothing(figsize) ? Dict{Symbol,Any}() : Dict{Symbol,Any}(:size => figsize)

    # Extract metadata about the glacier
    lon = results.lon
    lat = results.lat
    x = results.x
    y = results.y
    rgi_id = results.rgi_id
    Δx = results.Δx

    velocity_vars = [:V, :V_ref] # Velocity variables

    for var in intersect(velocity_vars, variables)  # Check only given vars for maximum
        if hasproperty(results, var)
            current_matrix = getfield(results, var)
            if !isnothing(current_matrix) && !isempty(current_matrix)
                if typeof(current_matrix) <: Vector
                    @assert length(current_matrix)>0 "Variable $(var) is an empty vector"
                    @assert (isnothing(timeIdx)) || (size(current_matrix,1)>=timeIdx) "The provided index=$(timeIdx) is greater than the size of the vector for $(var) which is $(size(current_matrix,1))"
                else
                end
            end
        end
    end

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
        ax = Axis(fig[1, ax_col], aspect=DataAspect())
        data = getfield(results, var)
        title, unit = get(title_mapping, string(var), (string(var), ""))

        if typeof(data) <: Vector
            @assert length(data)>0 "Variable $(var) is an empty vector"
            @assert (isnothing(timeIdx)) || (size(data,1)>=timeIdx) "The provided index=$(timeIdx) is greater than the size of the vector for $(var) which is $(size(data,1))"
            data = isnothing(timeIdx) ? data[end] : data[timeIdx]
            Vx = getfield(results, var == :V ? :Vx : :Vx_ref)
            Vy = getfield(results, var == :V ? :Vy : :Vy_ref)
            dataVx = isnothing(timeIdx) ? Vx[end] : Vx[timeIdx]
            dataVy = isnothing(timeIdx) ? Vy[end] : Vy[timeIdx]
        end

        X,Y=Sleipnir.meshgrid(x,y)

        positions = Point2f.(reshape(X,:), reshape(Y,:))
        directions = Vec2f.(dataVx, -dataVy)
        arrows2d!(ax, positions, directions; tiplength=tiplength, lengthscale=lengthscale)

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
    tspan::Tuple{F,F}=results.tspan,
    metrics::Vector{String}="difference",
    figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing,
) where {F<:AbstractFloat}
    # Check if more than one variable is passed
    @assert length(variables) == 1 "Only one variable can be passed to this function."

    figKwargs = isnothing(figsize) ? Dict{Symbol,Any}() : Dict{Symbol,Any}(:size => figsize)

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
    @assert typeof(data) == Vector{Matrix{Float64}} "Only temporal quantities can be used in this function."

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
            xlabel="Δ$variable_title ($(title_mapping[string(variables[1])][2]))",
            ylabel="Frequency",
            title="Histogram of $(title_mapping[string(variables[1])][1]) evolution ($rgi_id)"
        )
    elseif metrics == ["difference"]
        fig = Figure(; figKwargs...)
        ax_diff = Axis(fig[1, 1], title="$(title_mapping[string(variables[1])][1]) evolution ($rgi_id)",aspect=DataAspect())
    else
        figKwargs[:layout] = GridLayout(1, 3)
        fig = Figure(; figKwargs...)
        ax = Axis(
            fig[1, 3],
            xlabel="Δ$variable_title ($(title_mapping[string(variables[1])][2]))",
            ylabel="Frequency",
            title="Histogram of $(title_mapping[string(variables[1])][1])\n evolution"
        )
        ax_diff = Axis(fig[1, 1], title="$(title_mapping[string(variables[1])][1]) evolution",aspect=DataAspect())
        fig[0, :] = Label(fig, "$rgi_id")
    end

    # Plot based on the metric
    for metric in metrics
        if metric == "hist"
            hist!(ax, vec(data_diff[mask]), bins=50)
            ax.limits[] = (minimum(data_diff[mask]), maximum(data_diff[mask]), nothing, nothing)
        elseif metric == "difference"

            nx, ny = size(data_diff)

            # Calculate the symmetric color range
            max_abs_value = max(abs(minimum(data_diff[mask])), abs(maximum(data_diff[mask])))

            hm_diff = heatmap!(ax_diff, reverseForHeatmap(data_diff, x, y), colormap=:redsblues, colorrange=(-max_abs_value, max_abs_value))

            ax_diff.xlabel = "Longitude"
            ax_diff.ylabel = "Latitude"
            ax_diff.xticks=([round(nx/2)], ["$(round(lon;digits=6)) °"])
            ax_diff.yticks=([round(ny/2)], ["$(round(lat;digits=6)) °"])
            ax_diff.yticklabelrotation = π/2
            ax_diff.ylabelpadding = 5
            ax_diff.yticklabelalign = (:center, :bottom)

            # Width of the scale division in heatmap data units
            scale_width = 0.10*nx
            scale_number = round(Δx * scale_width / 1000; digits=1)#to km

            if metrics == ["difference"]
                textsize=1.2*scale_width
            else
                textsize=0.9*scale_width
            end

            # Position and draw the scale division rectangle
            poly!(ax_diff, Rect(nx -round(0.15*nx) , round(0.075*ny), scale_width, scale_width/10), color=:black)
            text!(ax_diff, "$scale_number km",
                position = (nx - round(0.15*nx)+scale_width/16, round(0.075*ny)+scale_width/10),
                fontsize=textsize)
            cb = Colorbar(fig[1, 2], hm_diff)
            Observables.connect!(cb.height, @lift CairoMakie.Fixed($(viewport(ax_diff.scene)).widths[2]))
            Label(fig[1, 2], "Δ$variable_title ($(title_mapping[string(variables[1])][2]))", fontsize=14, valign=:top, padding=(0, -25))

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
    metrics="median",
    tspan,
    threshold=0.5,
    figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing,
)
    # Check if more than one variable is passed
    @assert length(variables) == 1 "Only one variable can be passed to this function."

    figKwargs = isnothing(figsize) ? Dict{Symbol,Any}() : Dict{Symbol,Any}(:size => figsize)

    # Extract data for the variable
    data = getfield(results, variables[1])

    # Check the shape of the extracted data
    if typeof(data) ≠ Vector{Matrix{Float64}}
        error("Only temporal quantities can be used in this function.")
    end

    # Check for valid metrics
    valid_metrics = ["average", "median", "max", "std","min"]
    for metric in metrics
        if !(metric in valid_metrics)
            error("Invalid metric: $metric. Valid metrics are: $valid_metrics")
            return
        end
    end

    # Extract the rgi_id
    rgi_id = results.rgi_id

    # Create a time vector
    t = range(tspan[1], stop=tspan[2], length=length(getfield(results, variables[1])))

    # Create a single plot for all other metrics
    fig = Figure(; figKwargs...)
    ax = Axis(fig[1, 1], xlabel="Time (years)", ylabel="$(title_mapping[string(variables[1])][1]) ($(title_mapping[string(variables[1])][2]))", title="Metrics for $(title_mapping[string(variables[1])][1]) through time ($rgi_id)")

    # If "average" or "std" is in metrics, calculate them
    if "average" in metrics || "std" in metrics
        avg_vals = [mean(filter(x -> !(isnan(x)) && x >= threshold, matrix[:])) for matrix in data]
        std_vals = [std(filter(x -> !(isnan(x)) && x >= threshold, matrix[:])) for matrix in data]
    end

    # Calculate and plot metrics
    for metric in metrics
        if metric == "average"
            if "std" in metrics
                band!(ax, t, avg_vals .- std_vals, avg_vals .+ std_vals, alpha=0.1, label="Std Dev", color=:lightgray)
            end
            lines!(ax, t, avg_vals, label="Average")
        elseif metric == "median"
            median_vals = [median(filter(x -> !isnan(x) && x >= threshold, matrix[:])) for matrix in data]
            lines!(ax, t, median_vals, label="Median")
        elseif metric == "min"
            min_vals = [minimum(filter(x -> !isnan(x) && x >= threshold, matrix[:])) for matrix in data]
            lines!(ax, t, min_vals, linestyle=:dot, label="Min")
        elseif metric == "max"
            max_vals = [maximum(filter(x -> !isnan(x) && x >= threshold, matrix[:])) for matrix in data]
            lines!(ax, t, max_vals, linestyle=:dot, label="Max")
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
    figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing,
)

    # Determine pixel area
    area=results.Δx*results.Δy

    # Check if more than one variable is passed
    @assert length(variables) == 1 "Only one variable can be passed to this function."

    figKwargs = isnothing(figsize) ? Dict{Symbol,Any}() : Dict{Symbol,Any}(:size => figsize)

    # Extract the rgi_id
    rgi_id = results.rgi_id

    # Create a time vector
    t = range(tspan[1], stop=tspan[2], length=length(getfield(results, variables[1])))

    # Extract data for the variable
    data = getfield(results, variables[1])

    # Check the shape of the extracted data
    if typeof(data) ≠ Vector{Matrix{Float64}}
        error("Only temporal quantities can be used in this function.")
    end

    # Calculate integrated ice volume for each time step
    integrated_ice_volume = [sum(matrix) * area for matrix in data]/(10^9)  # Multiply by area and convert to km^3

    nDigitsYticks = max(2, Int(round(-log10(maximum(integrated_ice_volume) - minimum(integrated_ice_volume)))))

    # Plot the integrated ice volume as a function of time
    fig = Figure(; figKwargs...)
    ax = Axis(fig[1, 1]; ytickformat = x -> string.(round.(x; digits=nDigitsYticks)))
    lines!(ax, t, integrated_ice_volume, color=:blue)
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
    figsize::Union{Nothing, Tuple{Int64, Int64}}=nothing,
)
    # Check for exactly two variables
    @assert length(variables) == 2 "Exactly two variables are required for the scatter plot."

    figKwargs = isnothing(figsize) ? Dict{Symbol,Any}() : Dict{Symbol,Any}(:size => figsize)

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
        mask_H = getfield(results, obs_key) .> treshold[1] * maximum(getfield(results, obs_key))
        mask_H_2 = getfield(results, obs_key) .< (1 - treshold[2]) * maximum(getfield(results, obs_key))
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
    ax = Axis(fig[1, 1], xlabel = string(variables[1]), ylabel = string(variables[2]), title = "Scatter plot for RGI ID: " * rgi_id)

    scatter!(ax, vec(x_values), vec(y_values), markersize = 5, color = :blue, label = "Data")
    xmin, xmax = minimum(x_values), maximum(x_values)
    ymin, ymax = minimum(y_values), maximum(y_values)
    lines!(ax, [xmin, xmax], [xmin, xmax], linestyle = :dash, color = :red, label = "y = x")

    # Display metrics on the plot
    metrics_text = "RMSE: $(round(rmse, digits=2))\nR²: $(round(r_squared, digits=2))\nBias: $(round(bias, digits=2))"
    text!(ax, metrics_text, position = (xmax, ymax), align = (:right, :top), color = :black)

    return fig
end


"""
    plot_glacier(results::Results, plot_type::String, variables::Vector{Symbol}; kwargs...) -> Figure

Generate various types of plots for glacier data.

# Arguments
- `results::Results`: The results object containing the data to be plotted.
- `plot_type::String`: Type of plot to generate. Options are:
  * "heatmaps": Heatmaps for glacier variables like `:H`, `:H₀`, `:S`, `:B`, `:V`, `:Vx`, `:Vy`, `:V_ref`.
  * "evolution difference": Temporal difference metrics (between start and end) for a variable, with optional metrics like "hist" (histogram) and "difference".
  * "evolution statistics": Temporal statistical metrics for a variable, with optional metrics like "average", "median", "min", "max", and "std".
  * "integrated volume": Temporal evolution of the integrated ice volume for a variable.
  * "bias": Scatter plot to visualize the bias between two variables.
- `variables::Vector{Symbol}`: Variables to be plotted, e.g., `:H`.

# Optional Keyword Arguments
- `tspan`: A tuple representing the start and end time for the simulation.
- `metrics`: Metrics to visualize, e.g., `["average"]` for statistics, `["difference"]` for difference.
- `scale_text_size::Union{Nothing,Float64}`: Optional argument to scale the text size for heatmaps.
- `threshold::Vector{F}`: Threshold values for filtering data in bias plots.
- `figsize::Tuple{Int64, Int64}`: Size of the figure.

# Returns
- A `Figure` object containing the desired visualization.

# Notes
- Ensure the `variables` and `kwargs` match the requirements of the specified `plot_type`.
- The function routes requests to specific plotting functions based on `plot_type`.
"""
function plot_glacier(results::Results, plot_type::String, variables::Vector{Symbol}; kwargs...)

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

    if plot_type == "heatmaps"
        return plot_glacier_heatmaps(results, variables, title_mapping; kwargs...)
    elseif plot_type == "quivers"
        return plot_glacier_quivers(results, variables, title_mapping; kwargs...)
    elseif plot_type == "evolution difference"
        return plot_glacier_difference_evolution(results, variables, title_mapping; kwargs...)
    elseif plot_type == "evolution statistics"
        return plot_glacier_statistics_evolution(results, variables, title_mapping; kwargs...)
    elseif plot_type == "integrated volume"
        return plot_glacier_integrated_volume(results, variables, title_mapping; kwargs...)
    elseif plot_type == "bias"
        return plot_bias(results, variables; kwargs...)
    else
        error("Invalid plot_type: $plot_type")
    end
end

min_non_zero(M::Matrix{<: AbstractFloat}) = minimum(M[(!isnan).(M) .& (M .> 0)])

function plot_gridded_data(
    gridded_data::Union{Vector{Matrix{F}}, Matrix{F}},
    results::Results;
    scale_text_size::Union{Nothing,Float64}=nothing,
    timeIdx::Union{Nothing,Int64}=nothing,
    figsize::Union{Nothing, Tuple{Int64, Int64}} = nothing,
    plotContour::Bool=false,
    colormap = :cool,
    logPlot = false,
) where {F <: AbstractFloat}
    figKwargs = isnothing(figsize) ? Dict{Symbol,Any}() : Dict{Symbol,Any}(:size => figsize)

    # Extract metadata about the glacier
    lon = results.lon
    lat = results.lat
    x = results.x
    y = results.y
    rgi_id = results.rgi_id
    Δx = results.Δx
    mask = results.H[begin] .> 0.0
    nx, ny = size(results.H[begin])
    ctr = plotContour ? Contour.contour(collect(1:nx),1+ny.-collect(1:ny),mask, 0.5) : nothing

    @assert !isempty(gridded_data) "There is no data."
    if typeof(gridded_data) <: Vector
        @assert length(gridded_data)>0 "Data is an empty vector"
        @assert (isnothing(timeIdx)) || (size(gridded_data,1)>=timeIdx) "The provided index=$(timeIdx) is greater than the size of the vector which is $(size(gridded_data,1))"
        min_values = isnothing(timeIdx) ? min_non_zero(gridded_data[end]) : minimum(
            map(gridded_data) do M
                min_non_zero(M)
            end
        )
        max_values = maximum(replace(isnothing(timeIdx) ? gridded_data[end] : gridded_data[timeIdx], NaN => 0.))
    else
        min_values = min_non_zero(gridded_data)
        max_values = maximum(replace(gridded_data, NaN => 0.))
    end

    # Determine global minimum/maximum
    global_min = isempty(min_values) ? nothing : maximum(min_values)
    global_max = isempty(max_values) ? nothing : maximum(max_values)

    figKwargs[:layout] = GridLayout(2, 2)
    fig = Figure(; figKwargs...)

    ax_row = 1
    ax_col = 1
    ax = Axis(fig[ax_row, ax_col], aspect=DataAspect())
    data = deepcopy(gridded_data)
    # title, unit = get(title_mapping, string(var), (string(var), ""))

    if typeof(data) <: Vector
        @assert length(data)>0 "Data is an empty vector"
        @assert (isnothing(timeIdx)) || (size(data,1)>=timeIdx) "The provided index=$(timeIdx) is greater than the size of the vector which is $(size(data,1))"
        data = isnothing(timeIdx) ? data[end] : data[timeIdx]
    end

    nx, ny = size(data)

    mask = results.H[begin] .> 0.0
    data[.!mask] .= NaN

    hm = heatmap!(ax, reverseForHeatmap(data, x, y), colormap=colormap, colorrange=(logPlot ? global_min*0.85 : 0, global_max), colorscale=logPlot ? log10 : identity)
    cb = Colorbar(fig[ax_row, ax_col + 1], hm)
    Observables.connect!(cb.height, @lift CairoMakie.Fixed($(viewport(ax.scene)).widths[2]))
    # Label(fig[ax_row, ax_col + 1], "$var ($unit)", fontsize=14, valign=:top, padding=(0, -25))

    if plotContour
        for curve in ctr.lines
            xs = first.(curve.vertices)
            ys = last.(curve.vertices)
            lines!(ax, xs, ys, color=:black, linewidth=1)
        end
    end

    # ax.title = "$title"
    ax.xlabel = "Longitude"
    ax.ylabel = "Latitude"
    ax.xticks=([round(nx/2)], ["$(round(lon;digits=6)) °"])
    ax.yticks=([round(ny/2)], ["$(round(lat;digits=6)) °"])
    ax.yticklabelrotation = π/2
    ax.ylabelpadding = 5
    ax.yticklabelalign = (:center, :bottom)

    scale_width = 0.10*nx
    scale_number = round(Δx * scale_width / 1000; digits=1) # Convert to km
    textsize = isnothing(scale_text_size) ? 1.2*scale_width : scale_text_size

    poly!(ax, Rect(nx - round(0.15*nx), round(0.075*ny), scale_width, scale_width/10), color=:black)
    text!(ax, "$scale_number km", position=(nx - round(0.15*nx) + scale_width/16, round(0.075*ny) + scale_width/10), fontsize=textsize)

    fig[0, :] = Label(fig, "$rgi_id")
    resize_to_layout!(fig)
    return fig
end
