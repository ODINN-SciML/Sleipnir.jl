
function topography_and_inputs()
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]

    params = Parameters(
        simulation = SimulationParameters(
        use_velocities = false,
        use_glathida_data = true,
        multiprocessing = false,
        workers = 1,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        rgi_paths = rgi_paths,
        climate_data_source = :W5E5
    )
    )

    glacier = initialize_glaciers(rgi_ids, params)[1]
    simulation = Simulation(params, [glacier])

    # Test iTopoSlope and iTopoAspect: construction, names, and get_input output
    @testset "Topographic Inputs (iTopoSlope, iTopoAspect)" begin
        slope_input = iTopoSlope(window_m = 100.0)
        aspect_input = iTopoAspect(window_m = 100.0)

        @test default_name(slope_input) == :slope
        @test default_name(aspect_input) == :aspect

        slope_field = get_input(slope_input, simulation, 1, 2011.5)
        @test size(slope_field) == size(glacier.S)
        @test all(0.0 .≤ slope_field .≤ 90.0)

        aspect_field = get_input(aspect_input, simulation, 1, 2011.5)
        @test size(aspect_field) == size(glacier.S)
        @test all(0.0 .≤ aspect_field .< 360.0)
    end

    # Test the surface topography computation pipeline with synthetic surfaces
    @testset "Surface Topography Computation" begin
        S = Float64.(reshape(1:100, 10, 10))
        Δx, Δy = 50.0, 50.0

        slope, aspect = compute_surface_topography(S, Δx, Δy; window_m = 100.0)
        @test size(slope) == size(S)
        @test compute_surface_slope(S, Δx, Δy) ≈ slope
        @test compute_surface_aspect(S, Δx, Δy) ≈ aspect
        @test all(0.0 .≤ aspect .< 360.0)

        # Flat surface: slope must be zero everywhere
        @test all(compute_surface_topography(fill(100.0, 10, 10), Δx, Δy)[1] .≈ 0.0)

        # Glacier2D dispatch
        slope_g, _ = compute_surface_topography(glacier)
        @test size(slope_g) == size(glacier.S)
    end

    # Test _era5_fields_present for different climate step configurations
    @testset "ERA5 Fields Detection" begin
        cs_zero = ClimateStep(
            prcp = 0.1, temp = -5.0, gradient = -0.006,
            albedo = 0.0, slhf = 0.0, sshf = 0.0, ssrd = 0.0, str = 0.0,
            avg_temp = 0.0, avg_gradient = -0.006, ref_hgt = 2500.0
        )
        cs_era5 = ClimateStep(
            prcp = 0.1, temp = -5.0, gradient = -0.006,
            albedo = 0.5, slhf = 100.0, sshf = 50.0, ssrd = 200.0, str = -100.0,
            avg_temp = 0.0, avg_gradient = -0.006, ref_hgt = 2500.0
        )
        @test !_era5_fields_present(cs_zero)
        @test _era5_fields_present(cs_era5)

        # For Climate2Dstep, check detection with matrix fields
        cs2d_zero = Climate2Dstep(
            temp = fill(-5.0, 3, 3), PDD = fill(0.0, 3, 3),
            snow = fill(0.1, 3, 3), rain = fill(0.05, 3, 3),
            elevation_diff = fill(0.0, 3, 3), aspect = fill(180.0, 3, 3),
            albedo = fill(0.0, 3, 3), slhf = fill(0.0, 3, 3),
            slope = fill(15.0, 3, 3), sshf = fill(0.0, 3, 3),
            ssrd = fill(0.0, 3, 3), str = fill(0.0, 3, 3),
            gradient = -0.006, avg_gradient = -0.006,
            x = collect(1:3), y = collect(1:3), ref_hgt = 2500.0
        )
        cs2d_era5 = Climate2Dstep(
            temp = fill(-5.0, 3, 3), PDD = fill(0.0, 3, 3),
            snow = fill(0.1, 3, 3), rain = fill(0.05, 3, 3),
            elevation_diff = fill(0.0, 3, 3), aspect = fill(180.0, 3, 3),
            albedo = fill(0.2, 3, 3), slhf = fill(100.0, 3, 3),
            slope = fill(15.0, 3, 3), sshf = fill(50.0, 3, 3),
            ssrd = fill(200.0, 3, 3), str = fill(-100.0, 3, 3),
            gradient = -0.006, avg_gradient = -0.006,
            x = collect(1:3), y = collect(1:3), ref_hgt = 2500.0
        )
        @test !_era5_fields_present(cs2d_zero)
        @test _era5_fields_present(cs2d_era5)
    end

    # Test utility functions: slicing (including error path) and aggregation fallback
    @testset "Climate Slicing and Aggregation Utilities" begin
        # Use the glacier's raw climate which was loaded during initialization
        climate = glacier.climate.raw_climate

        climate_slice = _slice_climate_between_dates(climate, Date(2010, 6, 1), Date(2010, 8, 31))
        @test size(climate_slice, Ti) > 0

        @test_throws ArgumentError _slice_climate_between_dates(
            climate, Date(1900, 1, 1), Date(1900, 12, 31))

        # Layer exists → non-zero result
        @test _aggregate_raw_layer(climate_slice, :temp) != 0.0
        # Layer absent → returns zero without error
        @test _aggregate_raw_layer(climate_slice, :nonexistent) == 0.0
    end
end
