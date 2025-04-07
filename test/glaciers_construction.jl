

function glaciers2D_constructor(; save_refs::Bool = false, use_glathida_data::Bool = false)

    rgi_paths = get_rgi_paths()
    if use_glathida_data
        rgi_ids = ["RGI60-07.00042", "RGI60-07.00065"] # Use glaciers that have glathida data
        file_suffix = "w_glathida"
    else
        rgi_ids = ["RGI60-11.03638", "RGI60-11.01450"]
        file_suffix = "wo_glathida"
    end
    # Filter out glaciers that are not used to avoid having references that depend on all the glaciers processed in Gungnir
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    params = Parameters(simulation=SimulationParameters(velocities=false,
                                                        use_glathida_data=use_glathida_data,
                                                        working_dir=Sleipnir.root_dir,
                                                        test_mode=true,
                                                        rgi_paths=rgi_paths))

    glaciers = initialize_glaciers(rgi_ids, params)

    # Test prints
    println(glaciers)
    println(glaciers[1])
    println(glaciers[1].climate)
    println(glaciers[1].climate.climate_2D_step)

    if save_refs
        jldsave(joinpath(Sleipnir.root_dir, string("test/data/glaciers/glaciers2D_", file_suffix, ".jld2")); glaciers)
    end

    glaciers_ref = load(joinpath(Sleipnir.root_dir, string("test/data/glaciers/glaciers2D_", file_suffix, ".jld2")))["glaciers"]

    @test all(glaciers == glaciers_ref)


end

