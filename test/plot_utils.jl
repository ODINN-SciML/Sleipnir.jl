function glaciers2D_plots()
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    params = Parameters(
        simulation = SimulationParameters(
        use_velocities = false,
        use_glathida_data = false,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        multiprocessing = false,
        workers = 1,
        rgi_paths = rgi_paths
    )
    )
    glaciers = initialize_glaciers(rgi_ids, params)

    S = glaciers[1].S
    ifm = SimpleIceflowModel{Sleipnir.Float}(S)

    results = Results(
        glaciers[1],
        ifm;
        H = [abs.(randn(size(glaciers[1].H₀)...)), abs.(randn(size(glaciers[1].H₀)...))],
        Vx = [abs.(randn(size(glaciers[1].H₀)...)), abs.(randn(size(glaciers[1].H₀)...))],
        Vy = [abs.(randn(size(glaciers[1].H₀)...)), abs.(randn(size(glaciers[1].H₀)...))],
        Vx_ref = [
            abs.(randn(size(glaciers[1].H₀)...)), abs.(randn(size(glaciers[1].H₀)...))],
        Vy_ref = [
            abs.(randn(size(glaciers[1].H₀)...)), abs.(randn(size(glaciers[1].H₀)...))]
    )

    # Test execution
    @testset "plot_glacier tests" begin
        @testset "Heatmaps" begin
            plot_glacier(results, "heatmaps", [:H, :B]; plotContour = true)
        end

        @testset "Quivers" begin
            plot = plot_glacier(results, "quivers", [:V_ref, :V])
        end

        @testset "Statistics Evolution" begin
            plot_glacier(results, "evolution statistics", [:H], tspan = (2010.0, 2015.0),
                metrics = ["average", "std", "max", "median", "min"])
        end

        @testset "Difference Evolution" begin
            plot_glacier(results, "evolution difference", [:H], tspan = (2010.0, 2015.0),
                metrics = ["difference", "hist"], figsize = (800, 600))
        end

        @testset "Integrated Volume" begin
            plot_glacier(results, "integrated volume", [:H], tspan = (2010.0, 2015.0))
        end

        @testset "Bias" begin
            plot_glacier(results, "bias", [:B, :S])
        end

        @testset "Gridded data" begin
            plot_gridded_data(abs.(randn(size(glaciers[1].H₀)...)), results; logPlot = true)
        end

        @testset "Cumulative mass balance" begin
            # No MB history -> graceful nothing
            @test plot_cumulative_mb(results) === nothing
            # With MB history -> a Figure (raw and annually-averaged)
            mb_maps = [abs.(randn(size(glaciers[1].H₀)...)) for _ in 1:3]
            results_mb = Results(
                glaciers[1], ifm; H = results.H, MB = mb_maps, tspan = (2010.0, 2015.0))
            @test plot_cumulative_mb(results_mb) isa Figure
            @test plot_cumulative_mb(results_mb; annual_MB = true) isa Figure
        end

        @testset "Glacier DEM" begin
            @test plot_glacier_dem(results) isa Figure
            @test plot_glacier_dem(glaciers[1]) isa Figure
        end

        @testset "Observations" begin
            # Zero means no observation in the thickness and velocity data
            footprint = glaciers[1].H₀ .> 0
            H_obs = [footprint .* (50.0 .+ 100.0 .* rand(size(footprint)...)) for _ in 1:2]
            V_obs = [footprint .* (10.0 .+ 40.0 .* rand(size(footprint)...))]
            glacier_obs = Glacier2D(
                glaciers[1];
                thicknessData = ThicknessData([2010.0, 2015.0], H_obs),
                velocityData = SurfaceVelocityData(
                    vabs = V_obs,
                    date = [DateTime(2012, 7, 1)],
                    date1 = [DateTime(2012, 1, 1)],
                    date2 = [DateTime(2013, 1, 1)]
                )
            )

            @test plot_observations(glacier_obs; tspan = (2009.0, 2016.0)) isa Figure
            @test plot_observations(glacier_obs; max_points = 10) isa Figure
            # No thickness or velocity data in this glacier
            @test_throws ArgumentError plot_observations(glaciers[1])
            @test plot_observation_timeline([glacier_obs, glaciers[1]];
                tspan = (2008.0, 2017.0),
                glacier_tspans = Dict(glacier_obs.rgi_id => (2009.0, 2016.0))) isa Figure

            # Model against the thickness observations, at each date and averaged
            results_t = Results(glaciers[1], ifm; H = results.H, t = [2010.0, 2015.0])
            @test plot_glacier_vs_observations(results_t, glacier_obs, :H) isa Figure
            @test plot_glacier_vs_observations(
                results_t, glacier_obs, :H; aggregate = :mean) isa Figure
            # No simulated velocity in these results
            @test_throws ArgumentError plot_glacier_vs_observations(
                results_t, glacier_obs, :V)
            # Velocity with its reference, and a mask of the cells
            results_v = Results(glaciers[1], ifm; H = results.H, V = results.H,
                V_ref = [results.H[1]], date_Vref = [2012.5], t = [2010.0, 2015.0])
            @test plot_glacier_vs_observations(results_v, glacier_obs, :V) isa Figure
            plot_gridded_data(abs.(randn(size(glaciers[1].H₀)...)), results;
                mask = glaciers[1].H₀ .> 0)
            @test plot_field_histogram(results.H[1];
                references = ["reference" => 1.0], logScale = true) isa Figure
        end

        @testset "Save figure" begin
            fig = plot_glacier(results, "heatmaps", [:H])
            path = joinpath(mktempdir(), "heatmap.png")
            @test save_figure(fig, path) == path
            @test isfile(path)
        end

        @testset "Quivers with missing reference velocity" begin
            # No Vx_ref/Vy_ref -> the V_ref panel is skipped, not a crash
            results_noref = Results(
                glaciers[1], ifm; H = results.H, Vx = results.Vx, Vy = results.Vy)
            @test plot_glacier(results_noref, "quivers", [:V, :V_ref]) isa Figure
        end
    end
end

# Fake iceflow model struct to test the plotting functions
struct IceflowModelTest{F <: AbstractFloat} <: AbstractModel
    S::Matrix{F}
end
function make_thickness_video_test()
    rgi_ids = ["RGI60-11.03646"]
    rgi_paths = get_rgi_paths()
    working_dir = joinpath(Sleipnir.root_dir, "test/data")

    δt = 1/12
    params = Parameters(
        simulation = SimulationParameters(
        use_MB = true,
        use_iceflow = true,
        use_velocities = true,
        use_glathida_data = false,
        tspan = (2014.0, 2015.0),
        step_MB = δt,
        working_dir = working_dir,
        multiprocessing = false,
        workers = 1,
        rgi_paths = rgi_paths,
        ice_thickness_source = :Farinotti19
    ),
    )

    glaciers = initialize_glaciers(rgi_ids, params)

    nSteps = (params.simulation.tspan[2]-params.simulation.tspan[1])/δt
    timeSteps = params.simulation.tspan[1] .+ collect(0:nSteps) .* δt
    H = [rand(glaciers[1].nx, glaciers[1].ny) for t in timeSteps]

    ifm = IceflowModelTest{Float64}(glaciers[1].S)
    results = Results(glaciers[1], ifm; H = H)

    # The output format is inferred from the file extension.
    mp4Path = mktempdir() * ".mp4"
    plot_glacier_vid("thickness", results, glaciers[1], params.simulation.tspan,
        δt, mp4Path; baseTitle = "Bossons glacier")
    @test isfile(mp4Path) && filesize(mp4Path) > 0

    gifPath = mktempdir() * ".gif"
    plot_glacier_vid("thickness", results, glaciers[1], params.simulation.tspan,
        δt, gifPath; baseTitle = "Bossons glacier")
    @test isfile(gifPath)
    @test startswith(String(read(gifPath)[1:6]), "GIF8")
end
