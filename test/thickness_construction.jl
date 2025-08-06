
function thickness_construction()
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]

    params = Parameters(
        simulation=SimulationParameters(
            use_velocities=false,
            use_glathida_data=true,
            working_dir=Sleipnir.root_dir,
            test_mode=true,
            rgi_paths=rgi_paths
        )
    )

    glacier = initialize_glaciers(rgi_ids, params)[1]
    t = [2010.0, 2010.1]
    H = [randn(glacier.nx, glacier.ny), randn(glacier.nx, glacier.ny)]

    thickness_data = ThicknessData(t, H)
    JET.@test_opt target_modules=(Sleipnir,) ThicknessData(t, H)

    glacier_with_thickness_data = Glacier2D(glacier, thicknessData=thickness_data)
    JET.@test_opt target_modules=(Sleipnir,) Glacier2D(glacier, thicknessData=thickness_data)

    @test thickness_data==glacier_with_thickness_data.thicknessData
end
