

function glaciers2D_constructor(; save_refs::Bool = false)

    rgi_ids = ["RGI60-11.03638", "RGI60-11.01450"]

    params = Parameters(simulation=SimulationParameters(velocities=false,
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

    @show glaciers[1] ≈ glaciers_ref[1]
    @show glaciers[2] ≈ glaciers_ref[2]

    @show glaciers[1].rgi_id ≈ glaciers_ref[1].rgi_id
    @show glaciers[1].gdir ≈ glaciers_ref[1].gdir
    @show glaciers[1].climate ≈ glaciers_ref[1].climate
    @show glaciers[1].H₀ ≈ glaciers_ref[1].H₀
    @show glaciers[1].S ≈ glaciers_ref[1].S
    @show glaciers[1].B ≈ glaciers_ref[1].B
    @show glaciers[1].V ≈ glaciers_ref[1].V
    @show glaciers[1].A ≈ glaciers_ref[1].A
    @show glaciers[1].C ≈ glaciers_ref[1].C
    @show glaciers[1].n ≈ glaciers_ref[1].n
    @show glaciers[1].slope ≈ glaciers_ref[1].slope
    @show glaciers[1].dist_border ≈ glaciers_ref[1].dist_border
    @show glaciers[1].S_coords ≈ glaciers_ref[1].S_coords
    @show glaciers[1].Δx ≈ glaciers_ref[1].Δx
    @show glaciers[1].Δy ≈ glaciers_ref[1].Δy
    @show glaciers[1].nx ≈ glaciers_ref[1].nx
    @show glaciers[1].ny ≈ glaciers_ref[1].ny

    # @show glaciers[1].rgi_id
    # @show glaciers_ref[1].rgi_id
    # @show glaciers[1].gdir
    # @show glaciers_ref[1].gdir
    # @show glaciers[1].climate
    # @show glaciers_ref[1].climate
    

    
    @test all(glaciers .≈ glaciers_ref)

end