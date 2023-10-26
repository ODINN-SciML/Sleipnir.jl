function plot_glacier_heatmaps(results, variables, title_mapping)

    # Dictionary of variable-specific colormaps
    colormap_mapping = Dict(
        "H" => :YlGnBu,
        "S" => :terrain,
        "B" => :terrain,
        "V" => :viridis,
        "Vx" => :viridis,
        "Vy" => :viridis
    )


    # Extract the rgi_id 
    rgi_id = :rgi_id in fieldnames(typeof(results)) ? results.rgi_id : "none"

    # Number of variables to plot determines the number of subplots
    num_vars = length(variables)

    # Determine the grid layout based on the number of variables
    rows, cols = if num_vars == 1
        2, 2  # Adjusted for extra row
    elseif num_vars == 2
        3, 2  # Adjusted for extra row
    elseif num_vars in [3, 4]
        3, 4  # Adjusted for extra row
    else
        error("Unsupported number of variables.")
    end

    # Create a figure with GridLayout
    fig = Figure(layout=GridLayout(rows, cols))

    # Iterate over the variables and create subplots
    for (i, var) in enumerate(variables)
        ax_row = div(i - 1, 2) + 1
        ax_col = 2 * (rem(i - 1, 2)) + 1
        ax = Axis(fig[ax_row, ax_col], aspect=DataAspect())

        # Extract data from results for variable
        data = getfield(results, var)
        
        # If the data is a 3D array, take only the last matrix
        if typeof(data) == Vector{Matrix{Float64}}
            data = data[end]
        end

        # Mirror data around both X and Y axes
        data = reverse(reverse(data, dims=1),dims=2)

        # Fetch the colormap for the current variable from the mapping
        colormap = get(colormap_mapping, string(var), :cool)  # Default to :cool if not found

        # Plot the heatmap for the current variable with its designated colormap
        hm = heatmap!(ax, data, colormap=colormap)
        
        if num_vars==2
            Colorbar(fig[ax_row, ax_col + 1], hm)
        
        else    
            Colorbar(fig[ax_row, ax_col + 1], hm)
        end

        # Set title, labels, and other attributes for current variable
        ax.title = "$(get(title_mapping, string(var), (string(var), ""))[1]) ($(get(title_mapping, string(var), (string(var), ""))[2]))" # Fetch title from mapping, default to var if not found
        ax.xlabel = "X coordinate"
        ax.ylabel = "Y coordinate"
    end
    
    fig[0, :] = Label(fig, "$rgi_id")

    return fig
end


function plot_glacier_difference_evolution(results, variables, title_mapping; tspan, metrics)
        
        # Check if more than one variable is passed
        if length(variables) > 1
            error("Only one variable can be passed to this function.")
        end
    
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
    
        # Check the shape of the extracted data
        if typeof(data) ≠ Vector{Matrix{Float64}}
            error("Only temporal quantities can be used in this function.")
        end
        
        # Extract the rgi_id 
        rgi_id = :rgi_id in fieldnames(typeof(results)) ? results.rgi_id : "none"
        
        # Print plot information
        variable_title = get(title_mapping, variables[1], variables[1])
              
    
        # Create a time vector
        t = range(tspan[1], stop=tspan[2], length=length(getfield(results, variables[1])))
    
        matrix_size = size(data[1])
        diff_width = 1.5 * matrix_size[1]
        diff_height = 1.5 * matrix_size[2]
        data_diff=data[end] - data[1]
     
         # Determine whether to create a single plot or a subplot
         if metrics == ["hist"]
             fig = Figure()
             ax = Axis(fig[1, 1], xlabel="Δ$variable_title ($(title_mapping[string(variables[1])][2]))", ylabel="Frequency", title="Histogram of $(title_mapping[string(variables[1])][1]) Evolution")
         elseif metrics == ["difference"]
             fig = Figure()
             ax_diff = Axis(fig[1, 1], title="$(title_mapping[string(variables[1])][1]) Evolution")
         else
             fig = Figure(layout=GridLayout(1, 4)) 
             ax = Axis(fig[1, 3:4], xlabel="Δ$variable_title ($(title_mapping[string(variables[1])][2]))", ylabel="Frequency", title="Histogram of $(title_mapping[string(variables[1])][1]) Evolution")
             ax_diff = Axis(fig[1, 1], title="$(title_mapping[string(variables[1])][1]) Evolution",width=diff_width,height=diff_height)
         end
    
        # Plot based on the metric
        for metric in metrics
            if metric == "hist"
                hist!(ax, vec(data[end]-data[1]), bins=50)
                ax.limits[] = (minimum(data_diff), maximum(data_diff), nothing, nothing)
            elseif metric == "difference"
                # Mirror data around x and y
                data_diff = reverse(reverse(data_diff, dims=1), dims=2)
                
                hm_diff = heatmap!(ax_diff, data_diff, colormap=:redsblues, color=:auto, halign=:right)
                ax_diff.xlabel = "X coordinate"
                ax_diff.ylabel = "Y coordinate"
                Colorbar(fig[1, 2], hm_diff)
              
                
            end
        end
        
        fig[0, :] = Label(fig, "$rgi_id")

        fig  # Return the main figure
    
end


function plot_glacier_statistics_evolution(results, variables, title_mapping; tspan, metrics, threshold=0.5)
    

    # Check if more than one variable is passed
    if length(variables) > 1
        error("Only one variable can be passed to this function.")
    end

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
    rgi_id = :rgi_id in fieldnames(typeof(results)) ? results.rgi_id : "none"
    
    # Print plot information
    variable_title = get(title_mapping, variables[1], variables[1])
    

    # Create a time vector
    t = range(tspan[1], stop=tspan[2], length=length(getfield(results, variables[1])))

    # Create a single plot for all other metrics
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Time (years)", ylabel="$(title_mapping[string(variables[1])][1]) ($(title_mapping[string(variables[1])][2]))", title="Metrics for $(title_mapping[string(variables[1])][1]) through Time ($rgi_id)")

    # If "average" or "std" is in metrics, calculate them
    if "average" in metrics || "std" in metrics
        avg_vals = [mean(filter(x -> !(isnan(x)) && x >= threshold, matrix[:])) for matrix in data]
        std_vals = [std(filter(x -> !(isnan(x)) && x >= threshold, matrix[:])) for matrix in data]
    end

    # Calculate and plot metrics
    for metric in metrics
        if metric == "average"
            if "std" in metrics
                band!(ax, t, avg_vals .- std_vals, avg_vals .+ std_vals, fillalpha=0.1, label="Std Dev", color=:lightgray)
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
    
    

    fig  # Return the main figure
end


function plot_glacier_integrated_volume(results, variables, title_mapping; tspan, area=1.0)
   

    # Check if more than one variable is passed
    if length(variables) > 1
        error("Only one variable can be passed to this function.")
    end
    
    # Extract the rgi_id 
    rgi_id = :rgi_id in fieldnames(typeof(results)) ? results.rgi_id : "none"

   
    # Print plot information
    variable_title = get(title_mapping, variables[1], variables[1])
        
    # Create a time vector
    t = range(tspan[1], stop=tspan[2], length=length(getfield(results, variables[1])))

    # Extract data for the variable
    data = getfield(results, variables[1])

    # Check the shape of the extracted data
    if typeof(data) ≠ Vector{Matrix{Float64}}
        error("Only temporal quantities can be used in this function.")
    end

    # Calculate integrated ice volume for each time step
    integrated_ice_volume = [sum(matrix) * area for matrix in data]  # Multiply by area

    # Plot the integrated ice volume as a function of time
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, t, integrated_ice_volume, color=:blue)
    ax.xlabel = "Time (years)"
    ax.ylabel = "Integrated Ice Volume (m³)   "
    ax.title = "Evolution of Integrated Ice Volume ($rgi_id)"

    return fig  # Return the main figure with the plot
end

"""
    plot_glacier(results::T, plot_type::String, variables::Vector{Symbol}; kwargs...) -> Figure

Generate various types of plots for glacier data.

# Arguments
- `results`: A custom type containing the results of a glacier simulation.
- `plot_type`: Type of plot to generate. Options are:
  * "heatmaps": Heatmaps for glacier variables like `:H`, `:S`, `:B`, `:V`, `:Vx`, and `:Vy`.
  * "difference": Temporal difference metrics (between start and end) for a variable, with optional metrics like "hist" (histogram) and "difference".
  * "statistics": Temporal statistical metrics for a variable, with optional metrics like "average", "median", "min", "max", and "std".
  * "integrated_volume": Temporal evolution of the integrated ice volume for a variable.
- `variables`: Variables to be plotted, e.g., `:H`.

# Optional Keyword Arguments
- `tspan`: A tuple representing the start and end time for the simulation.
- `metrics`: Metrics to visualize, e.g., `["average"]` for statistics, `["difference"]` for difference.
- `area`: Area (in square meters) of each pixel for integrated volume calculation. Default is `1.0`.

# Returns
- A `Figure` object containing the desired visualization.

# Notes
- Ensure the `variables` and `kwargs` match the requirements of the specified `plot_type`.
- The function routes requests to specific plotting functions based on `plot_type`.
"""
function plot_glacier(results::T, plot_type::String, variables::Vector{Symbol}; kwargs...) where T
    title_mapping = Dict(
        "H" => ("Ice Thickness", "m", :YlGnBu),
        "S" => ("Surface Topography", "m", :terrain),
        "B" => ("Bed Topography", "m", :terrain),
        "V" => ("Ice Surface Velocity", "m/s", :viridis),
        "Vx" => ("Ice Surface Velocity (X-direction)", "m/s", :viridis),
        "Vy" => ("Ice Surface Velocity (Y-direction)", "m/s", :viridis)
    )

    if plot_type == "heatmaps"
        return plot_glacier_heatmaps(results, variables, title_mapping)
    elseif plot_type == "evolution_difference"
        return plot_glacier_difference_evolution(results, variables, title_mapping; kwargs...)
    elseif plot_type == "evolution_statistics"
        return plot_glacier_statistics_evolution(results, variables, title_mapping; kwargs...)
    elseif plot_type == "integrated_volume"
        return plot_glacier_integrated_volume(results, variables, title_mapping; kwargs...)
    else
        error("Invalid plot_type: $plot_type")
    end
end

