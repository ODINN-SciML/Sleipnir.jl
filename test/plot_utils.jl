function glaciers2D_plots()

    # Load glacier data
    @load (joinpath(@__DIR__,"data/glaciers/glaciers2D_test_temporal.jld2")) results

    # Test execution
    @testset "plot_glacier tests" begin
        @testset "Heatmaps" begin
            try
                plot_glacier(results[2], "heatmaps", [:H,:B])
                @test true
            catch
                @test false
            end
        end

        @testset "Statistics Evolution" begin
            try
                plot_glacier(results[2], "evolution_statistics", [:H], tspan=(2010.0,2015.0), metrics=["average","std","max","median","min"])
                @test true
            catch
                @test false
            end
        end

        @testset "Difference Evolution" begin
            try
                plot_glacier(results[2], "evolution_difference", [:H], tspan=(2010.0,2015.0), metrics=["difference","hist"])
                @test true
            catch
                @test false
            end
        end

        @testset "Integrated Volume" begin
            try
                plot_glacier(results[2], "integrated_volume", [:H], tspan=(2010.0,2015.0))
                @test true
            catch
                @test false
            end
        end
    end
end






