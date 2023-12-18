
export oggm_config

struct OGGMparameters <: AbstractParameters
    working_dir::String
    paths::Union{PyDict, Nothing}
    params::Union{PyDict, Nothing}
    multiprocessing::Bool
    workers::Int64
    ice_thickness_source::String
    base_url::String
end

"""
    OGGMparameters(;
        working_dir::String = joinpath(homedir(), "OGGM/OGGM_data"),
        paths::Union{PyDict, Nothing} = nothing,
        paths::Union{PyDict, Nothing} = nothing,
        multiprocessing::Bool = false,
        workers::Int64 = 1,
        base_url::String = "https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.6/L1-L2_files/elev_bands/"
        )
Initializes OGGM and it configures its parameters.
Keyword arguments
=================
    - `working_dir`: Working directory were all the files will be stored.
    - `paths`: Dictionary for OGGM-related paths.
    - `params`: Dictionary for OGGM-related parameters.
    - `multiprocessing`: Determines if multiprocessing is used for OGGM.
    - `workers`: How many workers are to be used for OGGM multiprocessing.
    - `ice_thickness_source`: Source for the ice thickness dataset. Either `Millan22` of `Farinotti19`.
    - `base_url`: Base URL to download all OGGM data.
"""
function OGGMparameters(;
            working_dir::String = joinpath(homedir(), "OGGM/OGGM_data"),
            paths::Union{PyDict, Nothing} = nothing,
            params::Union{PyDict, Nothing} = nothing,
            multiprocessing::Bool = false,
            workers::Int64 = 1,
            ice_thickness_source::String = "Farinotti19",
            base_url::String = "https://cluster.klima.uni-bremen.de/~oggm/gdirs/oggm_v1.6/L1-L2_files/elev_bands/",
            test = false
            )

    @assert ((ice_thickness_source == "Millan22") || (ice_thickness_source == "Farinotti19")) "Wrong ice thickness source! Should be either `Millan22` or `Farinotti19`."

    # Build the OGGM parameters and configuration
    OGGM_parameters = OGGMparameters(working_dir, paths, params,
                                    multiprocessing, workers, 
                                    ice_thickness_source,
                                    base_url)

    return OGGM_parameters
end

Base.:(==)(a::OGGMparameters, b::OGGMparameters) = a.working_dir == b.working_dir && a.paths == b.paths && a.params == b.params && 
                                      a.multiprocessing == b.multiprocessing && a.workers == b.workers && a.ice_thickness_source == b.ice_thickness_source &&
                                      a.base_url == b.base_url

"""
    oggm_config()

Configures the basic paths and parameters for OGGM.
"""
function oggm_config(working_dir=joinpath(homedir(), "OGGM/OGGM_data"); oggm_processes=1)
    scope = @__MODULE__ # Capture current module to allow use from external packages (e.g. Huginn, Muninn and ODINN)
    @eval begin
    @everywhere begin
    @eval $scope begin
    cfg.initialize() # initialize OGGM configuration
    
    PATHS = PyDict(cfg."PATHS")  # OGGM PATHS
    PATHS["working_dir"] = $working_dir # Choose own custom path for the OGGM data
    PARAMS = PyDict(cfg."PARAMS")
    PARAMS["hydro_month_nh"]=1
    PARAMS["dl_verify"] = false
    PARAMS["continue_on_error"] = true # avoid stopping when a task fails for a glacier (e.g. lack of data)

    # Multiprocessing 
    multiprocessing = $oggm_processes > 1 ? true : false
    PARAMS["use_multiprocessing"] = multiprocessing # Let's use multiprocessing for OGGM
    if multiprocessing
        PARAMS["mp_processes"] = $oggm_processes
    end



    end # @eval Sleipnir
    end # @everywhere
    end # @eval

end