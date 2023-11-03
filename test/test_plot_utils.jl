function glaciers2D_plots()

    # Load glacier data
    @load (joinpath(@__DIR__,"data/glaciers/glaciers2D_test_temporal.jld2")) results
    
    function compare_fig_to_image(fig, reference_img_path)
        # Temporary save figure
        temp_filename = tempname() * ".png"
        save(temp_filename, fig)
    
        # Step 2: Read Images
        img1 = load(temp_filename)
        img2 = load(joinpath(@__DIR__,reference_img_path))
    
        # Perform a pixel-by-pixel comparison
        are_identical = all(img1 .== img2)
    
        # Clean up the temporary file
        rm(temp_filename)
    
        return are_identical
    end
    
    
    # Test execution
    @testset "plot_glacier tests" begin
        @testset "Heatmaps" begin
            fig = plot_glacier(results[2], "heatmaps", [:H,:B])
            @test isa(fig, Figure)
            @test compare_fig_to_image(fig, "data/figures/reference_plot_heatmaps.png")
        end

        @testset "Statistics Evolution" begin
            fig = plot_glacier(results[2], "evolution_statistics", [:H], tspan=(2010.0,2015.0), metrics=["average","std","max","median","min"])
            @test isa(fig, Figure)
            @test compare_fig_to_image(fig, "data/figures/reference_plot_evolution_statistics.png")
        end

        @testset "Difference Evolution" begin
            fig = plot_glacier(results[2], "evolution_difference", [:H], tspan=(2010.0,2015.0), metrics=["difference","hist"])
            @test isa(fig, Figure)
            @test compare_fig_to_image(fig, "data/figures/reference_plot_evolution_difference.png")
        end

        @testset "Integrated Volume" begin
            fig = plot_glacier(results[2], "integrated_volume", [:H], tspan=(2010.0,2015.0))
            @test isa(fig, Figure)
            @test compare_fig_to_image(fig, "data/figures/reference_plot_integrated_volume.png")
        end
    end
end





