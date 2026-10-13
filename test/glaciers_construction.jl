
function glaciers2D_simple_constructor()
    glacier = Glacier2D()
    println(glacier)
end

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

    params = Parameters(
        simulation = SimulationParameters(
        use_velocities = false,
        use_glathida_data = use_glathida_data,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        rgi_paths = rgi_paths
    )
    )
    JET.@test_opt target_modules=(Sleipnir,) Parameters(
        simulation = SimulationParameters(
        use_velocities = false,
        use_glathida_data = use_glathida_data,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        rgi_paths = rgi_paths
    )
    )

    glaciers = initialize_glaciers(rgi_ids, params)
    # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)
    # JET.@test_opt broken=true target_modules=(Sleipnir,) initialize_glaciers(
    #     rgi_ids, params)

    # Test prints
    println(glaciers)
    println(glaciers[1])
    println(glaciers[1].climate)
    println(glaciers[1].climate.climate_2D_step)

    if use_glathida_data
        # The dated campaigns are attached and listed with the other observations
        @test all(!isnothing(glacier.thicknessData) for glacier in glaciers)
        @test count(row -> row.kind == :H, observations(glaciers)) ==
              sum(length(glacier.thicknessData.t) for glacier in glaciers)
        @test issorted([row.t for row in observations(glaciers[1])])
    end
    # The time window of a glacier is the global one unless it has its own
    id = glaciers[1].rgi_id
    windows = params.simulation.glacier_tspans
    @test glacier_tspan(params, glaciers[1]) == params.simulation.tspan
    windows[id] = (2010.5, 2014.0)
    @test glacier_tspan(params, glaciers[1]) == (2010.5, 2014.0)
    @test_logs check_glacier_tspans(params, glaciers)
    windows[id] = (2009.0, 2014.0) # starts before tspan
    @test_throws ArgumentError check_glacier_tspans(params, glaciers)
    delete!(windows, id)
    @test_logs check_glacier_tspans(params, glaciers)
    # A window for an RGI ID with no glacier only warns
    windows["RGI60-00.00000"] = params.simulation.tspan
    @test_logs (:warn,) check_glacier_tspans(params, glaciers)
    delete!(windows, "RGI60-00.00000")
    # Windows given as integers are converted, and are shown with the parameters
    sp = SimulationParameters(working_dir = params.simulation.working_dir,
        glacier_tspans = Dict(id => (2010, 2012)))
    @test sp.glacier_tspans[id] === (Sleipnir.Float(2010), Sleipnir.Float(2012))
    @test occursin("own window for 1 glacier(s)", sprint(show, sp))

    # Survey dates at most 30 days apart are one campaign, also as a chain
    dates = Date.(["2009-03-01", "2009-03-20", "2009-04-15", "2011-04-01", "2009-03-20"])
    @test Sleipnir._glathida_campaign_index(dates, 30) == [1, 1, 1, 2, 1]

    if use_glathida_data
        # A start close to the first survey is moved onto it, a survey far away is ignored
        t_first = minimum(glaciers[1].thicknessData.t)
        windows[id] = (t_first - 0.02, t_first + 1.0)
        @test_logs (:warn,) align_start!(params, glaciers[1:1]; within = 0.05, apply = false)
        @test windows[id] == (t_first - 0.02, t_first + 1.0)
        @test_logs (:info,) align_start!(params, glaciers[1:1]; within = 0.05)
        @test windows[id] == (t_first, t_first + 1.0)
        windows[id] = (t_first - 5.0, t_first + 1.0)
        @test_logs align_start!(params, glaciers[1:1]; within = 0.05)
        @test windows[id] == (t_first - 5.0, t_first + 1.0)
        delete!(windows, id)
    end

    if save_refs
        jldsave(
            joinpath(Sleipnir.root_dir,
                string("test/data/glaciers/glaciers2D_", file_suffix, ".jld2"));
            glaciers)
    end

    glaciers_ref = load(joinpath(
        Sleipnir.root_dir, string("test/data/glaciers/glaciers2D_", file_suffix, ".jld2")))["glaciers"]

    if !all(glaciers == glaciers_ref)
        println("Variables glaciers and glaciers_ref are different")
        for i in 1:size(glaciers, 1)
            println("Glacier n°$i")
            println("Glacier fields identical: ",
                Sleipnir.diffToDict(glaciers[i], glaciers_ref[i]))
            if !(glaciers[i].climate == glaciers_ref[i].climate)
                println("Climate fields identical: = ",
                    Sleipnir.diffToDict(glaciers[i].climate, glaciers_ref[i].climate))
                if !(glaciers[i].climate.climate_2D_step ==
                     glaciers_ref[i].climate.climate_2D_step)
                    println("Climate 2D step fields identical: = ",
                        Sleipnir.diffToDict(glaciers[i].climate.climate_2D_step,
                            glaciers_ref[i].climate.climate_2D_step))
                end
            end
        end
    end
    @test all(glaciers .≈ glaciers_ref) # Test approximate equality
    @test all(glaciers .== glaciers_ref) # Test exact equality
end

function glaciers2D_thickness_source(ice_thickness_source::Symbol)
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-11.03638", "RGI60-11.01450"]
    # Filter out glaciers that are not used to avoid having references that depend on all the glaciers processed in Gungnir
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    params = Parameters(
        simulation = SimulationParameters(
        use_velocities = true,
        ice_thickness_source = ice_thickness_source,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        rgi_paths = rgi_paths
    )
    )
    JET.@test_opt target_modules=(Sleipnir,) Parameters(
        simulation = SimulationParameters(
        use_velocities = true,
        ice_thickness_source = ice_thickness_source,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        rgi_paths = rgi_paths
    )
    )

    glaciers = initialize_glaciers(rgi_ids, params)
    # For the moment this is not type stable because of the readings (type of CSV files and RasterStack cannot be determined at compilation time)
    # JET.@test_opt broken=true target_modules=(Sleipnir,) initialize_glaciers(
    #     rgi_ids, params)

    # The Millan22 velocities point downhill: along each axis of the matrix, the velocity and
    # the surface slope have opposite signs (this fails with the wrong sign of vy)
    for glacier in glaciers
        ice = glacier.H₀[2:(end - 1), 2:(end - 1)] .> 10
        dSdi = glacier.S[3:end, 2:(end - 1)] .- glacier.S[1:(end - 2), 2:(end - 1)]
        dSdj = glacier.S[2:(end - 1), 3:end] .- glacier.S[2:(end - 1), 1:(end - 2)]
        @test mean((glacier.Vx[2:(end - 1), 2:(end - 1)] .* dSdi)[ice]) < 0
        @test mean((glacier.Vy[2:(end - 1), 2:(end - 1)] .* dSdj)[ice]) < 0
    end

    # With the Millan22 velocity product, the velocities are attached with their error (if the
    # error layers are in the data)
    params_V = Parameters(
        simulation = SimulationParameters(
        use_velocities = true,
        velocity_product = :Millan22,
        ice_thickness_source = ice_thickness_source,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        rgi_paths = rgi_paths
    )
    )
    velocityData = initialize_glaciers(rgi_ids, params_V)[1].velocityData
    @test !isnothing(velocityData)
    @test isnothing(velocityData.vabs_error) || all(isfinite, velocityData.vabs_error)

    # Also test to load the data directly because this is used in ODINN
    for rgi_id in rgi_ids
        if ice_thickness_source==:Farinotti19
            farinotti19_thickness(rgi_id, params)
        elseif ice_thickness_source==:Millan22
            millan22_thickness(rgi_id, params)
        end
    end
end
