

"""

A structure to hold simulation parameters for a simulation in ODINN.

    struct SimulationParameters{I <: Integer, F <: AbstractFloat, VM <: VelocityMapping} <: AbstractParameters

# Fields
- `use_MB::Bool`: Flag to indicate whether mass balance should be used.
- `use_iceflow::Bool`: Flag to indicate whether ice flow should be used.
- `plots::Bool`: Flag to indicate whether plots should be generated.
- `velocities::Bool`: Flag to indicate whether velocities should be calculated.
- `overwrite_climate::Bool`: Flag to indicate whether to overwrite climate data.
- `use_glathida_data::Bool`: Flag to indicate whether to use GLATHIDA data.
- `tspan::Tuple{F, F}`: Time span for the simulation.
- `step::F`: Time step for the simulation.
- `multiprocessing::Bool`: Flag to indicate whether multiprocessing should be used.
- `workers::I`: Number of workers for multiprocessing.
- `working_dir::String`: Directory for working files.
- `test_mode::Bool`: Flag to indicate whether to run in test mode.
- `rgi_paths::Dict{String, String}`: Dictionary of RGI paths.
- `ice_thickness_source::String`: Source of ice thickness data.
- `mapping::VM`: Mapping to use in order to grid the data from the coordinates of
    the velocity product datacube to the glacier grid.
"""
struct SimulationParameters{I <: Integer, F <: AbstractFloat, VM <: VelocityMapping} <: AbstractParameters
    use_MB::Bool
    use_iceflow::Bool
    plots::Bool
    velocities::Bool
    overwrite_climate::Bool
    use_glathida_data::Bool
    tspan::Tuple{F, F}
    step::F
    multiprocessing::Bool
    workers::I
    working_dir::String
    test_mode::Bool
    rgi_paths::Dict{String, String}
    ice_thickness_source::String
    mapping::VM
end



"""
Constructor for `SimulationParameters` type, including default values.

    SimulationParameters(;
        use_MB::Bool = true,
        use_iceflow::Bool = true,
        plots::Bool = true,
        velocities::Bool = true,
        overwrite_climate::Bool = false,
        use_glathida_data::Bool = false,
        tspan::Tuple{F, F} = (2010.0,2015.0),
        step::F = 1/12,
        multiprocessing::Bool = true,
        workers::I = 4,
        working_dir::String = "",
        test_mode::Bool = false,
        rgi_paths::Dict{String, String} = Dict{String, String}(),
        ice_thickness_source::String = "Farinotti19",
        mapping::VM = MeanDateVelocityMapping(),
    ) where {I <: Integer, F <: AbstractFloat, VM <: VelocityMapping}


# Keyword arguments
- `use_MB::Bool`: Whether to use mass balance (default: `true`).
- `use_iceflow::Bool`: Whether to use ice flow (default: `true`).
- `plots::Bool`: Whether to generate plots (default: `true`).
- `velocities::Bool`: Whether to calculate velocities (default: `true`).
- `overwrite_climate::Bool`: Whether to overwrite climate data (default: `false`).
- `use_glathida_data::Bool`: Whether to use GLATHIDA data (default: `false`).
- `float_type::DataType`: Data type for floating point numbers (default: `Float64`).
- `int_type::DataType`: Data type for integers (default: `Int64`).
- `tspan::Tuple{F, F}`: Time span for the simulation (default: `(2010.0, 2015.0)`).
- `step::F`: Time step for the simulation (default: `1/12`).
- `multiprocessing::Bool`: Whether to use multiprocessing (default: `true`).
- `workers::I`: Number of workers for multiprocessing (default: `4`).
- `working_dir::String`: Working directory for the simulation (default: `""`).
- `test_mode::Bool`: Whether to run in test mode (default: `false`).
- `rgi_paths::Dict{String, String}`: Dictionary of RGI paths (default: `Dict{String, String}()`).
- `ice_thickness_source::String`: Source of ice thickness data, either `"Millan22"` or `"Farinotti19"` (default: `"Farinotti19"`).
- `mapping::VM`: Mapping to use in order to grid the data from the coordinates of
    the velocity product datacube to the glacier grid.

# Returns
- `simulation_parameters`: A new `SimulationParameters` object.

# Throws
- `AssertionError`: If `ice_thickness_source` is not `"Millan22"` or `"Farinotti19"`.

# Notes
- If the global variable ODINN_OVERWRITE_MULTI is set to true, multiprocessing is
    disabled in any case. This is to fix the documentation generation as for the
    moment Literate.jl freezes when multiprocessing is enabled.
"""
function SimulationParameters(;
    use_MB::Bool = true,
    use_iceflow::Bool = true,
    plots::Bool = true,
    velocities::Bool = true,
    overwrite_climate::Bool = false,
    use_glathida_data::Bool = false,
    tspan::Tuple{F, F} = (2010.0,2015.0),
    step::F = 1/12,
    multiprocessing::Bool = true,
    workers::I = 4,
    working_dir::String = "",
    test_mode::Bool = false,
    rgi_paths::Dict{String, String} = Dict{String, String}(),
    ice_thickness_source::String = "Farinotti19",
    mapping::VM = MeanDateVelocityMapping(),
) where {I <: Integer, F <: AbstractFloat, VM <: VelocityMapping}

    @assert ((ice_thickness_source == "Millan22") || (ice_thickness_source == "Farinotti19")) "Wrong ice thickness source! Should be either `Millan22` or `Farinotti19`."

    # Literate.jl fails at generating the documentation when multiprocessing is enabled
    if lowercase(get(ENV, "ODINN_OVERWRITE_MULTI", "false")) == "true"
        workers = 1
        multiprocessing = false
    end

    simulation_parameters = SimulationParameters(use_MB, use_iceflow, plots, velocities,
                                                overwrite_climate, use_glathida_data,
                                                Sleipnir.Float.(tspan), Sleipnir.Float(step),
                                                multiprocessing, Sleipnir.Int(workers), working_dir,
                                                test_mode, rgi_paths, ice_thickness_source, mapping)

    if !ispath(working_dir)
        mkpath(joinpath(working_dir, "data"))
    end

    return simulation_parameters
end

Base.:(==)(a::SimulationParameters, b::SimulationParameters) = a.use_MB == b.use_MB && a.use_iceflow == b.use_iceflow && a.plots == b.plots &&
                                      a.velocities == b.velocities && a.overwrite_climate == b.overwrite_climate && a.use_glathida_data == b.use_glathida_data &&
                                      a.tspan == b.tspan && a.step == b.step && a.multiprocessing == b.multiprocessing &&
                                      a.workers == b.workers && a.working_dir == b.working_dir && a.test_mode == b.test_mode && a.rgi_paths == b.rgi_paths &&
                                      a.ice_thickness_source == b.ice_thickness_source && a.mapping == b.mapping
