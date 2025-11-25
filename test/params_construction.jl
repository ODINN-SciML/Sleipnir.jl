
function params_constructor_specified(; save_refs::Bool = false)
    rgi_id = "RGI60-11.03638"
    rgi_paths = get_rgi_paths()
    # Filter out glaciers that are not used to avoid having references that depend on all the glaciers processed in Gungnir
    rgi_paths = Dict(rgi_id => rgi_paths[rgi_id])

    physical_params = PhysicalParameters(
        ρ = 900.0,
        g = 9.81,
        ϵ = 1e-10,
        η₀ = 1.0,
        maxA = 8e-17,
        minA = 8.5e-20,
        maxTlaw = 1.0,
        minTlaw = -25.0,
        noise_A_magnitude = 5e-18
    )
    JET.@test_opt PhysicalParameters(
        ρ = 900.0,
        g = 9.81,
        ϵ = 1e-10,
        η₀ = 1.0,
        maxA = 8e-17,
        minA = 8.5e-20,
        maxTlaw = 1.0,
        minTlaw = -25.0,
        noise_A_magnitude = 5e-18
    )

    simulation_params = SimulationParameters(
        use_MB = true,
        use_iceflow = true,
        plots = false,
        use_velocities = false,
        overwrite_climate = false,
        use_glathida_data = false,
        tspan = (2010.0, 2015.0),
        multiprocessing = false,
        workers = 10,
        working_dir = "",
        rgi_paths = rgi_paths
    )
    JET.@test_opt SimulationParameters(
        use_MB = true,
        use_iceflow = true,
        plots = false,
        use_velocities = false,
        overwrite_climate = false,
        use_glathida_data = false,
        tspan = (2010.0, 2015.0),
        multiprocessing = false,
        workers = 10,
        working_dir = "",
        rgi_paths = rgi_paths
    )

    params = Parameters(
        physical = physical_params,
        simulation = simulation_params
    )
    JET.@test_opt target_modules=(Sleipnir,) Parameters( # Report only issues in Sleipnir
        physical = physical_params,
        simulation = simulation_params
    )

    if save_refs
        jldsave(
            joinpath(
                Sleipnir.root_dir, "test/data/params/simulation_params_specified.jld2");
            simulation_params)
        jldsave(
            joinpath(Sleipnir.root_dir, "test/data/params/physical_params_specified.jld2");
            physical_params)
        jldsave(
            joinpath(Sleipnir.root_dir, "test/data/params/params_specified.jld2"); params)
    end

    simulation_params_ref = load(joinpath(
        Sleipnir.root_dir, "test/data/params/simulation_params_specified.jld2"))["simulation_params"]
    physical_params_ref = load(joinpath(
        Sleipnir.root_dir, "test/data/params/physical_params_specified.jld2"))["physical_params"]
    params_ref = load(joinpath(
        Sleipnir.root_dir, "test/data/params/params_specified.jld2"))["params"]

    @test physical_params == physical_params_ref
    @test simulation_params == simulation_params_ref
    @test params == params_ref
end

function params_constructor_default(; save_refs::Bool = false)
    physical_params = PhysicalParameters()
    JET.@test_opt PhysicalParameters()

    simulation_params = SimulationParameters()
    JET.@test_opt SimulationParameters()

    params = Parameters(
        simulation = simulation_params,
        physical = physical_params
    )
    JET.@test_opt target_modules=(Sleipnir,) Parameters( # Report only issues in Sleipnir
        simulation = simulation_params,
        physical = physical_params
    )

    if save_refs
        jldsave(
            joinpath(Sleipnir.root_dir, "test/data/params/simulation_params_default.jld2");
            simulation_params)
        jldsave(
            joinpath(Sleipnir.root_dir, "test/data/params/physical_params_default.jld2");
            physical_params)
        jldsave(joinpath(Sleipnir.root_dir, "test/data/params/params_default.jld2"); params)
    end

    simulation_params_ref = load(joinpath(
        Sleipnir.root_dir, "test/data/params/simulation_params_default.jld2"))["simulation_params"]
    physical_params_ref = load(joinpath(
        Sleipnir.root_dir, "test/data/params/physical_params_default.jld2"))["physical_params"]
    params_ref = load(joinpath(Sleipnir.root_dir, "test/data/params/params_default.jld2"))["params"]

    @test physical_params == physical_params_ref
    @test simulation_params == simulation_params_ref
    @test params == params_ref
end
