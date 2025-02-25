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
    @load (joinpath(@__DIR__,"data/prediction/results2D_test.jld2")) results
    @load (joinpath(@__DIR__,"data/prediction/glaciers2D_test.jld2")) glaciers
    @load (joinpath(@__DIR__,"data/prediction/simuparams2D_test.jld2")) simulation

    tempPath = mktempdir()*".mp4"

    plot_glacier_vid("thickness", results[1], glaciers[1], simulation, tempPath; baseTitle="Bossons glacier")
end
