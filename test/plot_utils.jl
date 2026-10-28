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
            # The outline of the observed cells: one cell has 4 sides, two neighbours share one
            single = falses(5, 5)
            single[2, 3] = true
            @test length(Sleipnir._mask_outline(single)) == 2 * 4
            pair = copy(single)
            pair[3, 3] = true
            @test length(Sleipnir._mask_outline(pair)) == 2 * 6
            @test isempty(Sleipnir._mask_outline(falses(5, 5)))
            plot_gridded_data(abs.(randn(size(glaciers[1].H₀)...)), results;
                mask = glaciers[1].H₀ .> 0)
            @test plot_field_histogram(results.H[1];
                references = ["reference" => 1.0], logScale = true) isa Figure

            # Plots to inspect an inversion: the dated DEMs have a cell without value
            dem_a = copy(glaciers[1].S)
            dem_b = dem_a .- 5.0
            dem_b[findfirst(footprint)] = NaN
            glacier_i = Glacier2D(glacier_obs;
                surfaceData = SurfaceData(; t = [2010.5, 2014.5], S = [dem_a, dem_b],
                    σ = [zero(dem_a), zero(dem_b)]))
            results_i = Results(glaciers[1], ifm;
                H = [footprint .* 200.0, footprint .* 180.0],
                V = [footprint .* 30.0, footprint .* 28.0], t = [2010.0, 2015.0])
            @test plot_surface_change(results_i, glacier_i) isa Figure
            @test plot_bed_maps(results_i, glacier_i) isa Figure
            @test plot_fit_to_observations(results_i, glacier_i) isa Figure
            # The dated DEMs are observations too: they are in the 3D plot and in the timeline
            @test plot_observations(glacier_i; tspan = (2009.0, 2016.0)) isa Figure
            @test plot_observations(
                Glacier2D(glaciers[1]; surfaceData = glacier_i.surfaceData)) isa Figure
            @test plot_observation_timeline([glacier_i, glaciers[1]];
                tspan = (2008.0, 2017.0)) isa Figure
            # Through the router, which needs the glacier
            @test plot_glacier(results_i, "bed", Symbol[]; glacier = glacier_i,
                title = "Bed") isa Figure
            @test_throws ArgumentError plot_glacier(results_i, "fit", Symbol[])
            # A symbol is in camel case and is the string with spaces; an unknown type lists the valid ones
            @test plot_glacier(results_i, :surfaceChange, Symbol[]; glacier = glacier_i) isa Figure
            @test plot_glacier(results_i, "surface change", Symbol[]; glacier = glacier_i) isa Figure
            @test plot_glacier(results_i, :heatmaps, [:H]) isa Figure
            @test plot_glacier(results_i, :integratedVolume, [:H]; tspan = (2010.0, 2015.0)) isa Figure
            @test_throws ArgumentError plot_glacier(results_i, "not a plot", [:H])
            @test_throws ArgumentError plot_glacier(results_i, :surface_change, Symbol[]; glacier = glacier_i)
            # A glacier without dated DEMs cannot have a surface change
            @test_throws ArgumentError plot_surface_change(results_i, glaciers[1])
            # A bed of the model that is not the one of the results
            @test plot_bed_maps(results_i, glacier_i; B = glaciers[1].B .- 10.0) isa Figure

            # The figures show dates, never decimal years
            @test Sleipnir._date_label(2009.0) == "1 Jan 2009"
            @test Sleipnir._date_label(2000.5) == "2 Jul 2000"
            @test Sleipnir._date_label(DateTime(2012, 7, 1)) == "1 Jul 2012"
            # The colour range is the full range of the data: nothing is clipped
            panel = MapPanel([-1.0 5.0; 0.0 NaN], "x"; symmetric = true)
            @test Sleipnir._panel_limits(panel) == (-5.0, 5.0)
            @test Sleipnir._panel_limits(MapPanel([2.0 7.0; 3.0 NaN], "x")) == (2.0, 7.0)
            @test Sleipnir._panel_limits(MapPanel(fill(NaN, 2, 2), "x")) == (0.0, 1.0)
            # The relief of a plane is zero, also at the margin of the glacier (where it is not drawn)
            plane = [2.0 * i + 0.5 * j for i in 1:8, j in 1:8]
            window = falses(8, 8)
            window[2:7, 2:7] .= true
            relief_plane = Sleipnir._cell_relief(plane, window)
            @test all(isnan, relief_plane[.!window])
            @test all(isnan, relief_plane[[2, 7], :])
            @test all(abs.(relief_plane[3:6, 3:6]) .< 1e-12)
            # A bump of one cell has a relief
            bump = zeros(8, 8)
            bump[4, 4] = 9.0
            @test Sleipnir._cell_relief(bump, window)[4, 4] ≈ 8.0
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
