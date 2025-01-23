

function glaciers2D_constructor(; save_refs::Bool = false)

    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-11.03638", "RGI60-11.01450"]

    params = Parameters(simulation=SimulationParameters(velocities=false,
                                                        use_glathida_data=false,
                                                        working_dir=Sleipnir.root_dir,
                                                        test_mode=true,
                                                        rgi_paths=rgi_paths))

    glaciers = initialize_glaciers(rgi_ids, params; test=true)

    if save_refs
        jldsave(joinpath(Sleipnir.root_dir, "test/data/glaciers/glaciers2D.jld2"); glaciers)
    end

    glaciers_ref = load(joinpath(Sleipnir.root_dir,"test/data/glaciers/glaciers2D.jld2"))["glaciers"]

    @test all(glaciers == glaciers_ref)


end

