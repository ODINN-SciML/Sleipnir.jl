

struct SimulationParameters{F <: AbstractFloat} <: AbstractParameters
    use_MB::Bool
    use_iceflow::Bool
    plots::Bool
    velocities::Bool 
    overwrite_climate::Bool
    use_glathida_data::Bool
    float_type::DataType
    int_type::DataType
    tspan::Tuple{F, F}
    step::F
    multiprocessing::Bool
    workers::Int 
    working_dir::String
end


"""
    SimulationParameters(;
                        use_MB::Bool = true,
                        use_iceflow::Bool = true,
                        plots::Bool = true,
                        velocities::Bool = true,
                        overwrite_climate::Bool = false,
                        use_glathida_data::Bool = false,
                        float_type::DataType = Float64,
                        int_type::DataType = Int64,
                        tspan::Tuple{F, F} = (2010.0,2015.0),
                        multiprocessing::Bool = true,
                        workers::Int = 4
        )
Initialize the parameters for a simulation.
Keyword arguments
=================
    - `use_MB`: Determines if surface mass balance should be used.
    - `plots`: Determines if plots should be made.
    - `overwrite_climate`: Determines if climate data should be overwritten
    - 'use_glathida_data': Determines if data from the Glathida data set should be used
"""
function SimulationParameters(;
            use_MB::Bool = true,
            use_iceflow::Bool = true,
            plots::Bool = true,
            velocities::Bool = true,
            overwrite_climate::Bool = false,
            use_glathida_data::Bool = false,
            float_type::DataType = Float64,
            int_type::DataType = Int64,
            tspan::Tuple{F, F} = (2010.0,2015.0),
            step::F = 1/12,
            multiprocessing::Bool = true,
            workers::Int = 4,
            working_dir::String = ""
            ) where {F <: AbstractFloat}

    simulation_parameters = SimulationParameters(use_MB, use_iceflow, plots, velocities,
                                                overwrite_climate, use_glathida_data,
                                                float_type, int_type,
                                                tspan, step, multiprocessing, workers, working_dir)

    if !ispath(working_dir)
        mkpath(joinpath(working_dir, "data"))
    end

    return simulation_parameters
end

Base.:(==)(a::SimulationParameters, b::SimulationParameters) = a.use_MB == b.use_MB && a.use_iceflow == b.use_iceflow && a.plots == b.plots && 
                                      a.velocities == b.velocities && a.overwrite_climate == b.overwrite_climate && a.use_glathida_data == b.use_glathida_data &&
                                      a.float_type == b.float_type && a.int_type == b.int_type &&
                                      a.tspan == b.tspan && a.step == b.step && a.multiprocessing == b.multiprocessing &&
                                      a.workers == b.workers && a.working_dir == b.working_dir
