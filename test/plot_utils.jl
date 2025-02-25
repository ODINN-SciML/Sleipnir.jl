function glaciers2D_plots()

    # Load glacier data
    @load (joinpath(@__DIR__,"data/glaciers/glaciers2D_test_temporal.jld2")) results

    # Test execution
    @testset "plot_glacier tests" begin
        @testset "Heatmaps" begin
            try
                plot_glacier(results[1], "heatmaps", [:H,:B])
                @test true
            catch
                @test false
            end
        end

        @testset "Statistics Evolution" begin
            try
                plot_glacier(results[1], "evolution statistics", [:H], tspan=(2010.0,2015.0), metrics=["average","std","max","median","min"])
                @test true
            catch
                @test false
            end
        end

        @testset "Difference Evolution" begin
            try
                plot_glacier(results[1], "evolution difference", [:H], tspan=(2010.0,2015.0), metrics=["difference","hist"])
                @test true
            catch
                @test false
            end
        end

        @testset "Integrated Volume" begin
            try
                plot_glacier(results[1], "integrated volume", [:H], tspan=(2010.0,2015.0))
                @test true
            catch
                @test false
            end
        end

        @testset "Bias" begin
            try
                plot_glacier(results[1], "bias", [:B,:S])
                @test true
            catch
                @test false
            end
        end
    end
end


function make_thickness_video_test()
    rgi_ids = ["RGI60-11.03646"]
    rgi_paths = get_rgi_paths()
    working_dir = joinpath(Sleipnir.root_dir, "test/data")

    params = Parameters(
        simulation = SimulationParameters(
            use_MB = true,
            use_iceflow = true,
            velocities = true,
            use_glathida_data = false,
            tspan = (2014.0, 2015.0),
            working_dir = working_dir,
            multiprocessing = true,
            workers = 1,
            rgi_paths = rgi_paths,
            ice_thickness_source = "Farinotti19",
        ),
    )

    glaciers = initialize_glaciers(rgi_ids, params)

    nSteps = (params.simulation.tspan[2]-params.simulation.tspan[1])/params.simulation.step
    timeSteps = params.simulation.tspan[1].+collect(0:nSteps).*params.simulation.step
    H = [rand(glaciers[1].nx, glaciers[1].ny) for t in timeSteps]

    tempPath = mktempdir()*".mp4"

    plot_glacier_vid("thickness", H, glaciers[1], params.simulation, tempPath; baseTitle="Bossons glacier")
end
