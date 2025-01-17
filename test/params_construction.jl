
function params_constructor_specified(; save_refs::Bool = false)

    rgi_paths = get_rgi_paths()
    rgi_paths = Dict(k => string(v) for (k,v) in pairs(rgi_paths)) # Convert Dict{String, Any} to Dict{String, String}

    physical_params = PhysicalParameters(ρ = 900.0,
                                        g = 9.81,
                                        ϵ = 1e-3,
                                        η₀ = 1.0, 
                                        maxA = 8e-17,
                                        minA = 8.5e-20,
                                        maxTlaw = 1.0,
                                        minTlaw = -25.0,
                                        noise_A_magnitude = 5e-18)

    simulation_params = SimulationParameters(use_MB = true,
                                            use_iceflow = true,
                                            plots = false,
                                            velocities = false,
                                            overwrite_climate = false,
                                            use_glathida_data = false,
                                            float_type = Float64,
                                            int_type = Int64,
                                            tspan = (2010.0,2015.0),
                                            multiprocessing = false,
                                            workers = 10,
                                            working_dir = "",
                                            rgi_paths = rgi_paths)


    params = Parameters(physical=physical_params,
                        simulation=simulation_params)

    if save_refs
        jldsave(joinpath(Sleipnir.root_dir, "test/data/params/simulation_params_specified.jld2"); simulation_params)
        jldsave(joinpath(Sleipnir.root_dir, "test/data/params/physical_params_specified.jld2"); physical_params)
        jldsave(joinpath(Sleipnir.root_dir, "test/data/params/params_specified.jld2"); params)
    end

    simulation_params_ref = load(joinpath(Sleipnir.root_dir, "test/data/params/simulation_params_specified.jld2"))["simulation_params"]
    physical_params_ref = load(joinpath(Sleipnir.root_dir, "test/data/params/physical_params_specified.jld2"))["physical_params"]
    params_ref = load(joinpath(Sleipnir.root_dir, "test/data/params/params_specified.jld2"))["params"]

    @test physical_params == physical_params_ref
    @test simulation_params == simulation_params_ref
    @test params == params_ref

end

function params_constructor_default(; save_refs::Bool = false)

    physical_params = PhysicalParameters()

    simulation_params = SimulationParameters()

    params = Parameters(simulation=simulation_params,
                        physical=physical_params
                        )

    if save_refs
        jldsave(joinpath(Sleipnir.root_dir, "test/data/params/simulation_params_default.jld2"); simulation_params)
        jldsave(joinpath(Sleipnir.root_dir, "test/data/params/physical_params_default.jld2"); physical_params)
        jldsave(joinpath(Sleipnir.root_dir, "test/data/params/params_default.jld2"); params)
    end

    simulation_params_ref = load(joinpath(Sleipnir.root_dir, "test/data/params/simulation_params_default.jld2"))["simulation_params"]
    physical_params_ref = load(joinpath(Sleipnir.root_dir, "test/data/params/physical_params_default.jld2"))["physical_params"]
    params_ref = load(joinpath(Sleipnir.root_dir, "test/data/params/params_default.jld2"))["params"]

    @test physical_params == physical_params_ref
    @test simulation_params == simulation_params_ref
    @test params == params_ref

end