function glaciers2D_plots()

    # Load glacier data
    @load "glaciers2D_test_temporal.jld2" results
    
    function compare_fig_to_image(fig, reference_img_path)
        # Save the generated figure to a temporary PNG file
        tmp_file_path = tempname() * ".png"
        Makie.save(tmp_file_path, fig)
    
        # Load the saved image
        img_data = load(tmp_file_path)
        
        # Load the reference image
        reference_img = load(joinpath(@__DIR__,"reference_img_path"))
    
        # Optionally, you can delete the temporary file
        rm(tmp_file_path)
    
        return all(img_data .== reference_img)
    end
    
    
    # Test execution
    @testset "plot_glacier tests" begin
        @testset "Heatmaps" begin
            fig = plot_glacier(results[2], "heatmaps", [:H, :V, :B])
            display(fig)
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

run_tests()



