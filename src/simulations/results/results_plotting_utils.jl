export plot_glacier

function plot_glacier_heatmaps(results, variables, title_mapping)

    # Dictionary of variable-specific colormaps
    colormap_mapping = Dict(key => value[3] for (key, value) in title_mapping)

    # Extract the rgi_id 
    rgi_id = :rgi_id in fieldnames(typeof(results)) ? results.rgi_id : "none"

    #Extract longitude and latitude 
    lon = if hasproperty(results, :lon)
        results.lon
    elseif hasproperty(results.gdir, :cenlon)
        results.gdir.cenlon
    else
        nothing
    end
    
    lat = if hasproperty(results, :lat)
        results.lat
    elseif hasproperty(results.gdir, :cenlat)
        results.gdir.cenlat
    else
        nothing
    end
    

    Δx = results.Δx 
    
    ice_thickness_vars = [:H, :H₀, :H_glathida, :H_pred, :H_obs] # Symbol representation of the variable names
    max_values = []
    for var in ice_thickness_vars
        if hasproperty(results, var)
            push!(max_values, maximum(getfield(results, var)))
        end
    end
    global_max = maximum(max_values)
    

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

        ny, nx = size(data)
        
        # Fix alignment of matrix
        data = reverse(data', dims=2)
        
        # Fetch the colormap for the current variable from the mapping
        colormap = get(colormap_mapping, string(var), :cool)  # Default to :cool if not found

        # Plot the heatmap for the current variable with its designated colormap
        if var in ice_thickness_vars
            # Plot with a uniform color range for ice thickness variables
            hm = heatmap!(ax, data, colormap=colormap, colorrange=(0, global_max))
            Colorbar(fig[ax_row, ax_col + 1], hm)
        else
            # Plot other variables without the uniform color range
            hm = heatmap!(ax, data, colormap=colormap)
            Colorbar(fig[ax_row, ax_col + 1], hm)
        end

        
        
           
        
        

        # Set title, labels, and other attributes for current variable
        title, unit = get(title_mapping, string(var), (string(var), ""))
        ax.title = "$title ($unit)"
        ax.xlabel = "Longitude"
        ax.ylabel = "Latitude"
        ax.xticks=([round(nx/2)], ["$lon °"])
        ax.yticks=([round(ny/2)], ["$lat °"])
        ax.yticklabelrotation = π/2
        ax.ylabelpadding = 15
        ax.yticklabelalign = (:center, :bottom)


        # Width of the scale division in heatmap data units
        scale_width = 0.10*nx
        scale_number = round(Δx * scale_width / 1000; digits=1)#to km
        if num_vars == 1
            textsize=1.2*scale_width 
        elseif num_vars == 2
            textsize=0.9*scale_width 
        else
            textsize=0.5*scale_width
        
        end
        
        # Position and draw the scale division rectangle
        poly!(ax, Rect(nx -round(0.15*nx) , round(0.075*ny), scale_width, scale_width/10), color=:black)
        text!(ax, "$scale_number km", 
            position = (nx - round(0.15*nx)+scale_width/16, round(0.075*ny)+scale_width/10),
            fontsize=textsize)
        
    end
    
    fig[0, :] = Label(fig, "$rgi_id")
    resize_to_layout!(fig)
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
        
        #Extract longitude and latitude 
        lon = hasproperty(results, :lon) ? results.lon : "none"
        lat = hasproperty(results, :lat) ? results.lat : "none"
        
        #pixel width
        Δx = results.Δx
    
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
        diff_width = 1.0 * matrix_size[2]
        diff_height = 1.0 * matrix_size[1]
        data_diff=data[end] - data[1]
     
         # Determine whether to create a single plot or a subplot
         if metrics == ["hist"]
             fig = Figure()
             ax = Axis(fig[1, 1], xlabel="Δ$variable_title ($(title_mapping[string(variables[1])][2]))", ylabel="Frequency", title="Histogram of $(title_mapping[string(variables[1])][1]) Evolution")
         elseif metrics == ["difference"]
             fig = Figure()
             ax_diff = Axis(fig[1, 1], title="$(title_mapping[string(variables[1])][1]) Evolution",aspect=DataAspect())
         else
             fig = Figure(layout=GridLayout(1, 4)) 
             ax = Axis(fig[1, 3:4], xlabel="Δ$variable_title ($(title_mapping[string(variables[1])][2]))", ylabel="Frequency", title="Histogram of $(title_mapping[string(variables[1])][1]) Evolution",width=diff_width,height=diff_height)
             ax_diff = Axis(fig[1, 1], title="$(title_mapping[string(variables[1])][1]) Evolution",aspect=DataAspect())
         end
    
        # Plot based on the metric
        for metric in metrics
            if metric == "hist"
                hist!(ax, vec(data[end]-data[1]), bins=50)
                ax.limits[] = (minimum(data_diff), maximum(data_diff), nothing, nothing)
            elseif metric == "difference"
                
                
                ny, nx = size(data_diff)
                data_diff = reverse(data_diff',dims=2) # Fix alignment
                
                # Calculate the symmetric color range
                max_abs_value = max(abs(minimum(data_diff)), abs(maximum(data_diff)))

                    
                hm_diff = heatmap!(ax_diff, data_diff, colormap=:redsblues, halign=:right, colorrange=(-max_abs_value, max_abs_value))

                ax_diff.xlabel = "Longitude"
                ax_diff.ylabel = "Latitude"
                ax_diff.xticks=([round(nx/2)], ["$lon °"])
                ax_diff.yticks=([round(ny/2)], ["$lat °"])
                ax_diff.yticklabelrotation = π/2
                ax_diff.ylabelpadding = 15.0
                ax_diff.yticklabelalign = (:center, :bottom)

                
                # Width of the scale division in heatmap data units
                scale_width = 0.10*nx
                scale_number = round(Δx * scale_width / 1000; digits=1)#to km
                
                if metrics == ["difference"]
                    textsize=1.2*scale_width
                else 
                    textsize=0.5*scale_width
                end
                
                # Position and draw the scale division rectangle
                poly!(ax_diff, Rect(nx -round(0.15*nx) , round(0.075*ny), scale_width, scale_width/10), color=:black)
                text!(ax_diff, "$scale_number km", 
                    position = (nx - round(0.15*nx)+scale_width/16, round(0.075*ny)+scale_width/10),
                    fontsize=textsize)
                Colorbar(fig[1, 2], hm_diff)
                
              
                
            end
        end
        
        fig[0, :] = Label(fig, "$rgi_id")
        resize_to_layout!(fig)
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
    resize_to_layout!(fig)
    

    fig  # Return the main figure
end


function plot_glacier_integrated_volume(results, variables, title_mapping; tspan)
    
    # Determine pixel area
    area=results.Δx*results.Δy
    
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

    resize_to_layout!(fig)
    return fig  # Return the main figure with the plot
end

function plot_bias(data, keys; treshold = [0, 0])
    # Check for exactly two keys
    if length(keys) != 2
        error("Exactly two keys are required for the scatter plot.")
    end

     # Ensure treshold is an array of length 2
    if length(treshold) == 1
        treshold = [treshold[1], treshold[1]]
    end
    
    # Extract data
    rgi_id = data.rgi_id
    x_values = getfield(data,keys[1])
    y_values = getfield(data,keys[2])
    
    # Filter non-zero observations if necessary
    if :H_obs in keys || :V_obs in keys
        obs_key = :H_obs in keys ? :H_obs : :V_obs
        #non_zero_indices = findall(getfield(data,obs_key) .!= 0)
        mask_H = getfield(data,obs_key) .> treshold[1] * maximum(getfield(data,obs_key))
        mask_H_2 = getfield(data,obs_key) .< (1-treshold[2]) * maximum(getfield(data,obs_key))
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
    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1], xlabel = string(keys[1]), ylabel = string(keys[2]), title = "Scatter Plot for RGI ID: " * rgi_id)

    scatter!(ax, x_values, y_values, markersize = 5, color = :blue, label = "Data")
    xmin, xmax = minimum(x_values), maximum(x_values)
    ymin, ymax = minimum(y_values), maximum(y_values)
    lines!(ax, [xmin, xmax], [xmin, xmax], linestyle = :dash, color = :red, label = "y = x")

    # Display metrics on the plot
    metrics_text = "RMSE: $(round(rmse, digits=2))\nR²: $(round(r_squared, digits=2))\nBias: $(round(bias, digits=2))"
    text!(ax, metrics_text, position = (xmax, ymax), align = (:right, :top), color = :black)

    fig
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

# Returns
- A `Figure` object containing the desired visualization.

# Notes
- Ensure the `variables` and `kwargs` match the requirements of the specified `plot_type`.
- The function routes requests to specific plotting functions based on `plot_type`.
"""
function plot_glacier(results::T, plot_type::String, variables::Vector{Symbol}; kwargs...) where T
    
    title_mapping = Dict(
        "H" => ("Ice Thickness", "m", :YlGnBu),
        "H₀" => ("Ice Thickness", "m", :YlGnBu),
        "H_glathida" => ("Ice Thickness (GlaThiDa)", "m", :YlGnBu),
        "S" => ("Surface Topography", "m", :terrain),
        "B" => ("Bed Topography", "m", :terrain),
        "V" => ("Ice Surface Velocity", "m/s", :viridis),
        "Vx" => ("Ice Surface Velocity (X-direction)", "m/s", :viridis),
        "Vy" => ("Ice Surface Velocity (Y-direction)", "m/s", :viridis),
        "H_pred" => ("Predicted Ice Thickness", "m", :YlGnBu),
        "H_obs" => ("Observed Ice Thickness", "m", :YlGnBu),
        "H_diff" => ("Ice Thickness Difference", "m", :RdBu),
        "V_pred" => ("Predicted Ice Surface Velocity", "m/s", :viridis),
        "V_obs" => ("Observed Ice Surface Velocity", "m/s", :viridis),
        "V_diff" => ("Ice Surface Velocity Difference", "m/s", :RdBu)
    )

    if plot_type == "heatmaps"
        return plot_glacier_heatmaps(results, variables, title_mapping)
    elseif plot_type == "evolution_difference"
        return plot_glacier_difference_evolution(results, variables, title_mapping; kwargs...)
    elseif plot_type == "evolution_statistics"
        return plot_glacier_statistics_evolution(results, variables, title_mapping; kwargs...)
    elseif plot_type == "integrated_volume"
        return plot_glacier_integrated_volume(results, variables, title_mapping; kwargs...)
    elseif plot_type == "bias"  
        return plot_bias(results, variables; kwargs...)
    else
        error("Invalid plot_type: $plot_type")
    end
end
