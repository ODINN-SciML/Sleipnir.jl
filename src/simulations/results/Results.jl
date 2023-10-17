
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


function plot_results(results, variables;)
    
    # Extract the rgi_id 
    rgi_id = results.rgi_id
    
    # Number of variables to plot determines the number of subplots
    num_vars = length(variables)
    
    # Determine the grid layout based on the number of variables
    rows, cols = if num_vars == 1
        1, 2
    elseif num_vars == 2
        2, 2
    elseif num_vars in [3, 4]
        2, 4
    #elseif num_vars in [5, 6]  #more than 5 gives bad visibility
        #3, 4
    else
        error("Unsupported number of variables.")
    end
    
    # Adjust the rows for the title
    total_rows = rows + 1

    # Create a figure with GridLayout
    fig = Figure(layout=GridLayout(total_rows, cols)) 

    # Add title at the top
    supertitle = Label(fig[1, :], "$(rgi_id)", textsize=30,halign=:center)
    
    # Iterate over the variables and create subplots
    for (i, var) in enumerate(variables)
        ax_row = div(i - 1, 2) + 2
        ax_col = 2 * (rem(i - 1, 2)) + 1
        ax = Axis(fig[ax_row, ax_col], aspect=DataAspect())
        
        # Extract data from results for variable
        data = getfield(results, var)

        # If the data is a 3D array, take only the last matrix
        if typeof(data) == Vector{Matrix{Float64}}
            data = data[end]
        end
        
        # Plot the heatmap for the current variable
        hm = heatmap!(ax, data, color=:viridis)
        Colorbar(fig[ax_row, ax_col + 1], hm)
        
        # Set title, labels, and other attributes for current variable
        ax.title = "$(var)"
        ax.xlabel = "X coordinate"
        ax.ylabel = "Y coordinate"
    end
    
    return fig
end

