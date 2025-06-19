
"""
    surface_velocity_data()

Test the initialization of ice surface velocity structure w/ and w/o glacier gridding.
It creates fake `RasterStack` objects that mimic the structure of the true netCDF files.
"""
function surface_velocity_data()
    @testset "Fake interpolated datacube w/o glacier gridding" begin
        fakeRasterStack = Sleipnir.fake_interpolated_datacube()
        initialize_surfacevelocitydata(fakeRasterStack)
        JET.@test_opt broken=true target_modules=(Sleipnir,) initialize_surfacevelocitydata(fakeRasterStack) # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)
    end
    @testset "Fake multi datacube w/o glacier gridding" begin
        fakeRasterStack = Sleipnir.fake_multi_datacube()
        initialize_surfacevelocitydata(fakeRasterStack)
        JET.@test_opt broken=true target_modules=(Sleipnir,) initialize_surfacevelocitydata(fakeRasterStack) # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)
    end

    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-11.03646"]
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)
    params = Parameters(
        simulation=SimulationParameters(
            velocities=true,
            use_glathida_data=false,
            working_dir=Sleipnir.root_dir,
            test_mode=true,
            rgi_paths=rgi_paths
        )
    )
    JET.@test_opt target_modules=(Sleipnir,) Parameters(
        simulation=SimulationParameters(
            velocities=true,
            use_glathida_data=false,
            working_dir=Sleipnir.root_dir,
            test_mode=true,
            rgi_paths=rgi_paths
        )
    )
    glaciers = initialize_glaciers(rgi_ids, params)
    # JET.@test_opt broken=true target_modules=(Sleipnir,) initialize_glaciers(rgi_ids, params) # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)

    @testset "Fake interpolated datacube w glacier gridding" begin
        fakeRasterStack = Sleipnir.fake_interpolated_datacube()
        initialize_surfacevelocitydata(fakeRasterStack; glacier=glaciers[1])
        JET.@test_opt broken=true target_modules=(Sleipnir,) initialize_surfacevelocitydata(fakeRasterStack; glacier=glaciers[1]) # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)
    end
    @testset "Fake multi datacube w/ glacier gridding" begin
        fakeRasterStack = Sleipnir.fake_multi_datacube()
        initialize_surfacevelocitydata(fakeRasterStack; glacier=glaciers[1])
        JET.@test_opt broken=true target_modules=(Sleipnir,) initialize_surfacevelocitydata(fakeRasterStack; glacier=glaciers[1]) # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)
    end
end
