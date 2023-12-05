

function glaciers2D_constructor(; save_refs::Bool = true)

    rgi_ids = ["RGI60-11.03638", "RGI60-11.01450"]

    params = Parameters(simulation=SimulationParameters(velocities=false,
                                                        use_glathida_data=true,
                                                        working_dir=Sleipnir.root_dir),
                        OGGM=OGGMparameters(ice_thickness_source="Farinotti19"))

    glaciers = initialize_glaciers(rgi_ids, params; test=true)

    # Empty all PyCall stuff to avoid issues
    for glacier in glaciers
        glacier.gdir = nothing
        glacier.climate = nothing
        glacier.S_coords = nothing
    end

    if save_refs
        jldsave(joinpath(Sleipnir.root_dir, "test/data/glaciers/glaciers2D.jld2"); glaciers)
    end

    glaciers_ref = load(joinpath(Sleipnir.root_dir,"test/data/glaciers/glaciers2D.jld2"))["glaciers"]

    @test all(glaciers .≈ glaciers_ref)


end

