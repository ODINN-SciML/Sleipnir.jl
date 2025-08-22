
function helpers()
    @testset "glacierName" begin
        glacierName("RGI60-07.00042")=="Bodleybreen"
    end

    @testset "safe_approx" begin
        @test safe_approx(nothing, nothing)
        @test !safe_approx(nothing, 1.0)
        @test safe_approx(1, 1+1e-8)
    end
end

function simulation_utils()
    u = nothing
    integrator = nothing
    tstops = collect(2010.0:1/12:2011.0)
    @testset "stop_condition_tstops" begin
        @test Sleipnir.stop_condition_tstops(u, 2010.0, integrator, tstops)
        @test !Sleipnir.stop_condition_tstops(u, 2010.1, integrator, tstops)
    end
end

function results_utils()
    y = 2012
    m = 7
    d = 2 # This should be middle of the year
    dt = DateTime(y, m, d)
    @testset "datetime_to_floatyear" begin
        t = Sleipnir.datetime_to_floatyear(dt)
        @assert t==y+0.5 "t=$(t) but dt=$(dt)"
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
