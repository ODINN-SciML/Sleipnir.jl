include("iceflow_def.jl")

function results_default(; save_refs::Bool = false)
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]
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
    glaciers = initialize_glaciers(rgi_ids, params)

    S = glaciers[1].S
    ifm = SimpleIceflowModel{Sleipnir.Float}(S)

    results = Results(glaciers[1], ifm)

    if save_refs
        jldsave(joinpath(Sleipnir.root_dir, string("test/data/results/results.jld2")); results)
    end

    results_ref = load(joinpath(Sleipnir.root_dir, string("test/data/results/results.jld2")))["results"]

    @test all(results == results_ref)
end
