
function climate_downscale()
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
    period = period = partial_year(Day, t - step):Day(1):partial_year(Day, t)
    @testset "get_cumulative_climate" begin
        climate_step = get_cumulative_climate(glacier.climate.raw_climate)
    end
    @testset "get_cumulative_climate!" begin
        get_cumulative_climate!(glacier.climate, period)
    end
    @testset "downscale_2D_climate" begin
        climate_2D_step = downscale_2D_climate(glacier.climate.climate_step, glacier.S, glacier.Coords)
    end
    @testset "downscale_2D_climate!" begin
        downscale_2D_climate!(glacier)
    end
end

function dummy_climate()
    climate = Sleipnir.DummyClimate2D(longterm_temps = [-2.0])
end
