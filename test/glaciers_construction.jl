
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

    # Test the image of the glacier, shown only in terminals with colors.
    # The number of rows of the image depends on the size of the terminal, so we try several sizes.
    for glacier in glaciers,
        dsize in [(20, 80), (21, 80), (24, 100), (25, 100), (30, 120), (31, 120)]

        out_color = sprint(show, glacier; context = (:color => true, :displaysize => dsize))
        @test occursin("\e[", out_color)
        # The text must start at the beginning of a line, after the image
        @test occursin(r"(^|\n)Glacier2D ", out_color)
    end
    @test !occursin("\e[", sprint(show, glaciers[1]))

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

    # Also test to load the data directly because this is used in ODINN
    for rgi_id in rgi_ids
        if ice_thickness_source==:Farinotti19
            farinotti19_thickness(rgi_id, params)
        elseif ice_thickness_source==:Millan22
            millan22_thickness(rgi_id, params)
        end
    end
end
