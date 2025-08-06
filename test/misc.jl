
function helpers()
    @testset "glacierName" begin
        glacierName("RGI60-07.00042")=="Bodleybreen"
    end
end

function operations_glacier_mask()
    params = Parameters(
        simulation=SimulationParameters(
            test_mode=true,
            rgi_paths=get_rgi_paths(),
        )
    )
    glacier = initialize_glaciers(["RGI60-07.00042"], params)[1]
    H = glacier.H₀

    @testset "is_in_glacier" begin
        is_in_glacier(H, 3)
    end
end

function glacier_grid_downscaling()
    params = Parameters(
        simulation=SimulationParameters(
            test_mode=true,
            rgi_paths=get_rgi_paths(),
            gridScalingFactor=4,
        )
    )
    glacier = initialize_glaciers(["RGI60-07.00042"], params)[1]
    @test glacier.nx==22
    @test glacier.ny==30
    @test size(glacier.H₀)==(glacier.nx,glacier.ny)
end
