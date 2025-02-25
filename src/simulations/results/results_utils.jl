

"""
    store_results!(simulation::SIM, glacier_idx::I, solution) where {SIM <: Simulation, I <: Int}

Store the results of a simulation of a single glacier into a `Results`.
"""
function create_results(simulation::SIM, glacier_idx::I, solution, loss=nothing; light=false, batch_id::Union{Nothing, I}=nothing) where {SIM <: Simulation, I <: Integer}
    # The solution contains all the steps including the intermediate ones
    # This results in solution having multiple values for a given time step, we select the last one of each time step
    nSteps = (simulation.parameters.simulation.tspan[2]-simulation.parameters.simulation.tspan[1])/simulation.parameters.simulation.step
    timeSteps = simulation.parameters.simulation.tspan[1].+collect(0:nSteps).*simulation.parameters.simulation.step
    solStepIndices = [findlast(==(val), solution.t) for val in timeSteps]

    H = light ? [solution.u[begin],solution.u[end]] : solution.u[solStepIndices]
    # Simulations using Reverse Diff require an iceflow model per glacier
    if isnothing(batch_id)
        iceflow_model = simulation.model.iceflow
    else
        iceflow_model = simulation.model.iceflow[batch_id]
    end
    if !isnothing(simulation.model.machine_learning)
        θ = simulation.model.machine_learning.θ
    else
        θ = nothing
    end

    results = Results(simulation.glaciers[glacier_idx], iceflow_model;
                      H = H,
                      S = iceflow_model.S,
                      B = simulation.glaciers[glacier_idx].B,
                      V = iceflow_model.V,
                      Vx = iceflow_model.Vx,
                      Vy = iceflow_model.Vy,
                      Δx = simulation.glaciers[glacier_idx].Δx,
                      Δy = simulation.glaciers[glacier_idx].Δy,
                      lon = simulation.glaciers[glacier_idx].cenlon,
                      lat = simulation.glaciers[glacier_idx].cenlat,
                      nx = simulation.glaciers[glacier_idx].nx,
                      ny = simulation.glaciers[glacier_idx].ny,  
                      tspan = simulation.parameters.simulation.tspan,
                      θ = θ,
                      loss = loss
                    )

    return results
end

"""
    save_results_file!(results_list::Vector{Results{F}}, simulation::SIM; path::Union{String,Nothing}=nothing) where {F <: AbstractFloat, SIM <: Simulation}

Save simulation results which are provided as a list of `Results` into a `.jld2` file.
This function also overrides the `results`` attribute of `simulation`.
"""
function save_results_file!(results_list::Vector{Results{F, I}}, simulation::SIM; path::Union{String,Nothing}=nothing) where {F <: AbstractFloat, I <: Int, SIM <: Simulation}
    # Create path for simulation results
    if isnothing(path)
        predictions_path = joinpath(dirname(Base.current_project()), "data/results/predictions")
    else
        predictions_path = path
    end
    if !ispath(predictions_path)
        mkpath(predictions_path)
    end

    simulation.results = results_list

    tspan = simulation.parameters.simulation.tspan
    nglaciers = length(simulation.glaciers)
    jldsave(joinpath(predictions_path, "prediction_$(nglaciers)glaciers_$tspan.jld2"); simulation.results)
end
