# Dictionary of variable-specific titles
title_mapping = Dict(
    "H" => ("Ice Thickness", "m"),
    "S" => ("Surface Topography", "m"),
    "B" => ("Bed Topography", "m"),
    "V" => ("Ice Surface Velocity", "m/s"),
    "Vx" => ("Ice Surface Velocity (X-direction)", "m/s"),
    "Vy" => ("Ice Surface Velocity (Y-direction)", "m/s")
)


"""
    plot_glacier_heatmaps(results::T, variables::Vector{Symbol}) -> Figure

Plot heatmaps for various glacier variables.

# Arguments
- `results`: A custom type containing the results of a glacier simulation.
  * Must have fields corresponding to the variables to be plotted.
- `variables`: A vector of symbols representing the variables to be plotted.
  * Possible variables: `:H`, `:S`, `:B`, `:V`, `:Vx`, and `:Vy` or any other variable that has the same datastructure.

# Returns
- A `Figure` object containing the heatmaps for the specified glacier variables.

# Notes
- The function uses predefined colormaps for each variable or a default for variables that are not in the dictionary.
- Titles are provided based on a predefined mapping, data that is not in the dictionary just get their symbol as title.
"""
function plot_glacier_heatmaps(results::T, variables::Vector{Symbol})::Figure where T

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

    # Print plot information
    println("Glacier Visualization Summary")
    println("======================================")
    println("RGI ID of the glacier: $rgi_id")
    println("Quantities being plotted:")
    for var in variables
        # Use title_mapping for a descriptive name, default to var if not found
        println("- $(get(title_mapping, string(var), string(var)))")
    end

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
        data = reverse(reverse(data, dims=1), dims=2)

        # Fetch the colormap for the current variable from the mapping
        colormap = get(colormap_mapping, string(var), :cool)  # Default to :cool if not found

        # Debugging print
        println("Plotting variable: $(get(title_mapping, string(var), string(var))[1]) with colormap: $colormap")

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

    return fig
end

"""
    plot_glacier_difference(
        results::T, 
        variables::Vector{Symbol}; 
        tspan::Tuple{Float64, Float64}, 
        metrics::Vector{String} = ["difference"]
    ) -> Figure

Plot the difference metrics for a specific glacier variable over a given time span.

# Arguments
- `results`: A custom type containing the results of a glacier simulation.
  * Must have fields corresponding to the variables to be plotted.
- `variables`: A vector of symbols representing the variable to be plotted. Only one variable is allowed.
  * Possible variables: `:H` or any other data that has time evolution in it.
- `tspan`: A tuple representing the start and end time for the simulation.
- `metrics`: A vector of strings indicating the type of metrics to visualize. Possible values are: "hist" (histogram) and "difference". Defaults to `["difference"]`.

# Returns
- A `Figure` object containing the difference metrics for the specified glacier variable and or a histogram showing the distribution.

# Notes
- The function uses predefined titles for each variable based on a mapping any other variable just get their symbol as title.
- Only one variable can be plotted at a time with this function.
"""
function plot_glacier_difference(results::T, variables::Vector{Symbol}; tspan::Tuple{Float64, Float64}, metrics::Vector{String} = ["difference"]) ::Figure where T 
        
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
        println("Glacier Visualization Summary")
        println("======================================")
        println("RGI ID of the glacier: $rgi_id")
        println("Quantity being plotted: $(title_mapping[string(variables[1])][1]) $variable_title")
        
    
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
             fig = Figure(layout=GridLayout(1, 3)) 
             ax = Axis(fig[1, 3], xlabel="Δ$variable_title ($(title_mapping[string(variables[1])][2]))", ylabel="Frequency", title="Histogram of $(title_mapping[string(variables[1])][1]) Evolution")
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
    
        fig  # Return the main figure
    
end

"""
    plot_glacier_statistics(
        results::T, 
        variables::Vector{Symbol}; 
        tspan::Tuple{Float64, Float64}, 
        metrics::Vector{String} = ["average"], 
        threshold_ratio::Float64 = 0.001
    ) -> Figure

Plot statistical metrics for a specific glacier variable over a given time span.

# Arguments
- `results`: A custom type containing the results of a glacier simulation.
  * Must have fields corresponding to the variables to be plotted.
- `variables`: A vector of symbols representing the variable to be plotted. Only one variable is allowed.
  * Possible variables: `:H`, or any other data that has time evolution in it.
- `tspan`: A tuple representing the start and end time for the simulation.
- `metrics`: A vector of strings indicating the statistical metrics to visualize. 
  * Possible values: "average", "median", "min", "max", "std". Defaults to `["average"]`.
- `threshold_ratio`: A float specifying the threshold ratio for filtering data. Data below this threshold is ignored in the statistical calculations.

# Returns
- A `Figure` object containing the statistical metrics for the specified glacier variable.

# Notes
- The function uses predefined titles for each variable based on a mapping, any other variable just get their symbol as title.
- Only one variable can be plotted at a time with this function.
- The threshold for filtering is calculated based on the maximum value in the data and the specified threshold_ratio.
"""
function plot_glacier_statistics(results::T, variables::Vector{Symbol}; tspan::Tuple{Float64, Float64}, metrics::Vector{String} = ["average"], threshold_ratio::Float64 = 0.001)::Figure where T
    

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
    valid_metrics = ["average", "median", "min", "max", "std"]
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
    println("Glacier Visualization Summary")
    println("======================================")
    println("RGI ID of the glacier: $rgi_id")
    println("Quantity being plotted: $(title_mapping[string(variables[1])][1]) $variable_title")
    println("Threshold for filtering: $variable_title > $threshold_ratio • max_matrix_value")
    
    println("Statistics being plotted:")
    for met in metrics
        println("- $met")
    end

    # Create a time vector
    t = range(tspan[1], stop=tspan[2], length=length(getfield(results, variables[1])))

    # Create a single plot for all other metrics
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Time (years)", ylabel="$(title_mapping[string(variables[1])][1]) ($(title_mapping[string(variables[1])][2]))", title="Metrics for $(title_mapping[string(variables[1])][1]) through Time")

    # Calculate statistics, filtering out values below the threshold
    threshold = threshold_ratio * maximum(map(m -> maximum(m), data))
    threshold=round(threshold,digits=2)
    println("Threshold value for $(title_mapping[string(variables[1])][1]) when calculating statistics: $variable_title = $threshold")

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

"""
    plot_glacier_IV(
        results::T, 
        variables::Vector{Symbol}; 
        tspan::Tuple{Float64, Float64}, 
        area::Float64 = 1.0
    ) -> Figure

Plot the temporal evolution of the integrated ice volume for a specific glacier variable.

# Arguments
- `results`: A custom type containing the results of a glacier simulation.
  * Must have fields corresponding to the variables to be plotted.
- `variables`: A vector of symbols representing the variable to be integrated and plotted. Only one variable is allowed.
  * Possible variables: `:H` or any other data that has time evolution in it.
- `tspan`: A tuple representing the start and end time for the simulation.
- `area`: A float specifying the area (in square meters) of each pixel. This is used to calculate the integrated volume. Default is `1.0`.

# Returns
- A `Figure` object containing the plot showing the temporal evolution of integrated ice volume.

# Notes
- The function uses predefined titles for each variable based on a mapping, any other variable just get their symbol as title.
- Only one variable can be plotted at a time with this function.
- The integrated ice volume is calculated by summing up the values of the variable for each time step and multiplying by the provided pixel area.
"""
function plot_glacier_IV(results::T, variables::Vector{Symbol}; tspan::Tuple{Float64, Float64}, area::Float64 = 1.0)::Figure where T
   

    # Check if more than one variable is passed
    if length(variables) > 1
        error("Only one variable can be passed to this function.")
    end
    
    # Extract the rgi_id 
    rgi_id = :rgi_id in fieldnames(typeof(results)) ? results.rgi_id : "none"

   
    # Print plot information
    variable_title = get(title_mapping, variables[1], variables[1])
    println("Glacier Visualization Summary")
    println("======================================")
    println("RGI ID of the glacier: $rgi_id")
    println("Quantity being integrated: $(title_mapping[string(variables[1])][1]) $variable_title")
    if area==1.0
        println("Pixel Area is not specified")
    else 
        println("Pixel Area = $area m²")
    end
    
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
    ax.title = "Temporal Evolution of Integrated Ice Volume"

    return fig  # Return the main figure with the plot
end


