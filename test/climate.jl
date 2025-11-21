
function climate_downscale(; save_refs::Bool=false)
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]

    params = Parameters(
        simulation=SimulationParameters(
            use_velocities=false,
            use_glathida_data=true,
            working_dir=Sleipnir.root_dir,
            test_mode=true,
            rgi_paths=rgi_paths
        )
    )

    glacier = initialize_glaciers(rgi_ids, params)[1]

    step = 1/12
    t = 2011.0
    period = partial_year(Day, t - step):Day(1):partial_year(Day, t)

    # Perform climate downscaling
    climate_step = get_cumulative_climate(glacier.climate.raw_climate)
    get_cumulative_climate!(glacier.climate, t, step)
    climate_2D_step = downscale_2D_climate(glacier.climate.climate_step, glacier.S, glacier.Coords)
    downscale_2D_climate!(glacier)

    JET.@test_opt broken=true target_modules=(Sleipnir,) get_cumulative_climate(glacier.climate.raw_climate)
    JET.@test_opt broken=false target_modules=(Sleipnir,) get_cumulative_climate!(glacier.climate, t, step)
    JET.@test_opt broken=false target_modules=(Sleipnir,) downscale_2D_climate!(glacier) 

    if save_refs
        jldsave(joinpath(Sleipnir.root_dir, string("test/data/climate/climate_step.jld2")); climate_step)
        climate_step_period = glacier.climate.climate_step
        jldsave(joinpath(Sleipnir.root_dir, string("test/data/climate/climate_step_period.jld2")); climate_step_period)
        jldsave(joinpath(Sleipnir.root_dir, string("test/data/climate/climate_2D_step.jld2")); climate_2D_step)
    end

    climate_step_ref = load(joinpath(Sleipnir.root_dir, "test/data/climate/climate_step.jld2"))["climate_step"]
    climate_step_period_ref = load(joinpath(Sleipnir.root_dir, "test/data/climate/climate_step_period.jld2"))["climate_step_period"]
    climate_2D_step_ref = load(joinpath(Sleipnir.root_dir, "test/data/climate/climate_2D_step.jld2"))["climate_2D_step"]

    @test climate_step == climate_step_ref
    @test glacier.climate.climate_step == climate_step_period_ref
    @test climate_2D_step == climate_2D_step_ref
    @test glacier.climate.climate_2D_step == climate_2D_step_ref

end

function dummy_climate()
    climate = Sleipnir.DummyClimate2D(longterm_temps_scalar = [-2.0], longterm_temps_gridded = [ -2.0 -1.5; -1.0 -0.5])
end