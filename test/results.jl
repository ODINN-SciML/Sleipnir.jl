include("iceflow_def.jl")

function results_default(; save_refs::Bool = false, useDatacube::Bool = false)
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-11.03646"]
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    params = Parameters(
        simulation=SimulationParameters(
            velocities=false,
            use_glathida_data=false,
            working_dir=Sleipnir.root_dir,
            test_mode=true,
            rgi_paths=rgi_paths
        )
    )
    @inferred Parameters(
        simulation=SimulationParameters(
            velocities=false,
            use_glathida_data=false,
            working_dir=Sleipnir.root_dir,
            test_mode=true,
            rgi_paths=rgi_paths
        )
    )
    if useDatacube
        fakeRasterStack = fake_interpolated_datacube()
        glaciers = initialize_glaciers(rgi_ids, params; velocityDatacubes=Dict{String, RasterStack}(rgi_ids[1] => fakeRasterStack))
        @inferred initialize_glaciers(rgi_ids, params; velocityDatacubes=Dict{String, RasterStack}(rgi_ids[1] => fakeRasterStack))
        prefix = "_vel"
    else
        glaciers = initialize_glaciers(rgi_ids, params)
        @inferred initialize_glaciers(rgi_ids, params)
        prefix = ""
    end

    S = glaciers[1].S
    ifm = SimpleIceflowModel{Sleipnir.Float}(S)
    @inferred SimpleIceflowModel{Sleipnir.Float}(S)

    results = Results(glaciers[1], ifm)
    @inferred Results(glaciers[1], ifm)

    if save_refs
        jldsave(joinpath(Sleipnir.root_dir, string("test/data/results/results$(prefix).jld2")); results)
    end

    results_ref = load(joinpath(Sleipnir.root_dir, string("test/data/results/results$(prefix).jld2")))["results"]

    @test all(results == results_ref)
end
