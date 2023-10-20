
mutable struct Results{F <: AbstractFloat} 
    rgi_id::String
    H::Vector{Matrix{F}}
    S::Matrix{F}
    B::Matrix{F}
    V::Matrix{F}
    Vx::Matrix{F}
    Vy::Matrix{F}
end


function Results(glacier::AbstractGlacier, ifm::IF;
        rgi_id::String = glacier.rgi_id,
        H::Vector{Matrix{F}} = Matrix{F}([]),
        S::Matrix{F} = zeros(F, size(ifm.S)),
        B::Matrix{F} = zeros(F, size(ifm.B)),
        V::Matrix{F} = zeros(F, size(ifm.V)),
        Vx::Matrix{F} = zeros(F, size(ifm.Vx)),
        Vy::Matrix{F} = zeros(F, size(ifm.Vy))
            ) where {F <: AbstractFloat, IF <: AbstractModel}

    # Build the results struct based on input values
    results = Results(rgi_id, H, S, B,
                      V, Vx, Vy)

    return results
end


function plot_results(results, variables)
    ## Still need to add geographical coords
    # Dictionary of variable-specific colormaps
    colormap_mapping = Dict(
        "H" => :YlGnBu,
        "S" => :terrain,
        "B" => :terrain,
        "V" => :viridis,
        "Vx" => :viridis,
        "Vy" => :viridis
    )

    # Dictionary of variable-specific titles
    title_mapping = Dict(
        "H" => "Ice thickness",
        "S" => "Surface topography",
        "B" => "Bed topography",
        "V" => "Ice surface velocity",
        "Vx" => "Ice surface velocity (X-direction)",
        "Vy" => "Ice surface velocity (Y-direction)"
    )

    # Extract the rgi_id 
    rgi_id = results.rgi_id

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
        println("Plotting variable: $(get(title_mapping, string(var), string(var))) with colormap: $colormap")

        # Plot the heatmap for the current variable with its designated colormap
        hm = heatmap!(ax, data, colormap=colormap)
        
        if num_vars==2
            Colorbar(fig[ax_row, ax_col + 1], hm)
        
        else    
            Colorbar(fig[ax_row, ax_col + 1], hm)
        end

        # Set title, labels, and other attributes for current variable
        ax.title = get(title_mapping, string(var), string(var))  # Fetch title from mapping, default to var if not found
        ax.xlabel = "X coordinate"
        ax.ylabel = "Y coordinate"
    end

    return fig
end

function plot_results_temporal(results, variables; tspan::Tuple{Float64, Float64}, metrics::Vector{String} = ["average"])
    ## Still need to add geographical coords
    # Dictionary of variable-specific titles
    title_mapping = Dict(
        "H" => "Ice thickness",
        "S" => "Surface topography",
        "B" => "Bed topography",
        "V" => "Ice surface velocity",
        "Vx" => "Ice surface velocity (X-direction)",
        "Vy" => "Ice surface velocity (Y-direction)"
    )

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
    

    # Extract the rgi_id 
    rgi_id = results.rgi_id
    
    # Print plot information
    variable_title = get(title_mapping, variables[1], variables[1])
    println("Glacier Visualization Summary")
    println("======================================")
    println("RGI ID of the glacier: $rgi_id")
    println("Quantity being plotted: $variable_title")
    
    println("Metrics being plotted:")
    for met in metrics
        println("- $met")
    end

    # Check for valid metrics
    valid_metrics = ["average", "median", "min", "max", "std", "difference"]
    for metric in metrics
        if !(metric in valid_metrics)
            println("Invalid metric: ", metric)
            return
        end
    end

    # Create a time vector
    t = range(tspan[1], stop=tspan[2], length=length(getfield(results, variables[1])))

    matrix_size = size(data[1])
    diff_width = 1.5*matrix_size[1]
    diff_height = 1.5*matrix_size[2]
    other_plot_height = 2*diff_width / 2
    other_plot_width = 2*diff_width

    # Determine whether to create a single plot or a subplot
    if metrics == ["difference"]
        # Create a single plot for the "difference" metric
        fig = Figure()
        ax_diff = Axis(fig[1, 1], title="Difference End-Start")
    elseif all(m -> m != "difference", metrics)
        # Create a single plot for all other metrics
        fig = Figure()
        ax = Axis(fig[1, 1], xlabel="Years", title="Metrics for $(title_mapping[string(variables[1])]) through Time")
    else
        # Create a subplot for both types of metrics
        fig = Figure(layout=GridLayout(2, 1))
        ax = Axis(fig[1, 1], xlabel="Years", title="Metrics for $(title_mapping[string(variables[1])]) through Time", height=other_plot_height, width=other_plot_width)
        ax_diff = Axis(fig[2, 1], title="Difference End-Start", height=diff_height, width=diff_width)
    end

    # If "average" or "std" is in metrics, calculate them
    if "average" in metrics || "std" in metrics
        avg_vals = [mean(filter(x -> !(isnan(x) || x == 0), matrix[:])) for matrix in data]
        std_vals = [std(filter(x -> !(isnan(x) || x == 0), matrix[:])) for matrix in data]
    end

    # Calculate and plot metrics
    for metric in metrics
        if metric == "average"
            if "std" in metrics
                band!(ax, t, avg_vals .- std_vals, avg_vals .+ std_vals, fillalpha=0.1, label="Std Dev", color=:lightgray)
            end
            lines!(ax, t, avg_vals, label="Average")
        elseif metric == "median"
            median_vals = [median(filter(x -> !(isnan(x) || x == 0), matrix[:])) for matrix in data]
            lines!(ax, t, median_vals, label="Median")
        elseif metric == "min"
            min_vals = [minimum(filter(x -> !(isnan(x) || x == 0), matrix[:])) for matrix in data]
            lines!(ax, t, min_vals, linestyle=:dot, label="Min")
        elseif metric == "max"
            max_vals = [maximum(filter(x -> !(isnan(x) || x == 0), matrix[:])) for matrix in data]
            lines!(ax, t, max_vals, linestyle=:dot, label="Max")
        elseif metric == "difference"
            
            #mirror data around x and y
            data_diff=data[end] - data[1]
            data_diff = reverse(reverse(data_diff, dims=1), dims=2)
            
            hm_diff = heatmap!(ax_diff, data[end] - data[1], colormap=:redsblues, color=:auto,halign=:right)
            ax_diff.xlabel = "X coordinate"
            ax_diff.ylabel = "Y coordinate"
            if metrics == ["difference"]
                Colorbar(fig[1, 2], hm_diff)  # Place colorbar to the right when only "difference" is plotted
            else
                Colorbar(fig[2, 2], hm_diff,halign=:left)  # Place colorbar below when other metrics are also plotted
            end
        end
    end
   
    # Display the legend for non-"difference" metrics
    if metrics != ["difference"]
        Legend(fig[1, 2], ax)  # Create a legend next to the metrics plot based on the main axis
    end

    fig  # Return the main figure with both plots
end

function plot_results_IV(results, variables; tspan::Tuple{Float64, Float64}, area::Float64=1.0)      
    # Dictionary of variable-specific titles
    title_mapping = Dict(
        "H" => "Ice thickness",
        "S" => "Surface topography",
        "B" => "Bed topography",
        "V" => "Ice surface velocity",
        "Vx" => "Ice surface velocity (X-direction)",
        "Vy" => "Ice surface velocity (Y-direction)"
    )

    # Check if more than one variable is passed
    if length(variables) > 1
        error("Only one variable can be passed to this function.")
    end
    
    # Extract the rgi_id 
    rgi_id = results.rgi_id
    
    # Print plot information
    variable_title = get(title_mapping, variables[1], variables[1])
    println("Glacier Visualization Summary")
    println("======================================")
    println("RGI ID of the glacier: $rgi_id")
    println("Quantity being plotted: $variable_title")
    
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
    lines!(ax, t, integrated_ice_volume, label="Integrated Ice Volume", color=:blue)
    ax.xlabel = "Time"
    ax.ylabel = "Integrated Ice Volume"
    ax.title = "Temporal Evolution of Integrated Ice Volume"

    return fig  # Return the main figure with the plot
end