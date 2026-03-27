
function _write_synthetic_era5_monthly(path::String)
    mkpath(dirname(path))
    start_date = Date(2008, 1, 1)
    end_date = Date(2018, 12, 1)
    dates = collect(start_date:Month(1):end_date)
    ntime = length(dates)

    ds = NCDatasets.NCDataset(path, "c")
    try
        NCDatasets.defDim(ds, "time", ntime)
        vtime = NCDatasets.defVar(ds, "time", Float64, ("time",))
        vtime.attrib["units"] = "days since 2008-01-01 00:00:00"
        vtime.attrib["calendar"] = "proleptic_gregorian"
        vtime[:] = Float64.(Dates.value.(dates .- start_date))

        # Variables required by Sleipnir for both basic and ERA5-aware workflows.
        vtemp = NCDatasets.defVar(ds, "temp", Float32, ("time",))
        vprcp = NCDatasets.defVar(ds, "prcp", Float32, ("time",))
        vgrad = NCDatasets.defVar(ds, "gradient", Float32, ("time",))
        vfal = NCDatasets.defVar(ds, "fal", Float32, ("time",))
        vslhf = NCDatasets.defVar(ds, "slhf", Float32, ("time",))
        vsshf = NCDatasets.defVar(ds, "sshf", Float32, ("time",))
        vssrd = NCDatasets.defVar(ds, "ssrd", Float32, ("time",))
        vstr = NCDatasets.defVar(ds, "str", Float32, ("time",))

        month_idx = Float32.(1:ntime)
        annual_cycle = Float32.(sin.(2.0f0 * Float32(pi) .* month_idx ./ 12.0f0))
        shoulder_cycle = Float32.(cos.(2.0f0 * Float32(pi) .* month_idx ./ 6.0f0))

        # Synthetic but realistic monthly-scale climatology at glacier altitude.
        vtemp[:] = -10.0f0 .+ 11.0f0 .* annual_cycle
        vprcp[:] = max.(0.0f0, 0.06f0 .+ 0.03f0 .* shoulder_cycle)
        vgrad[:] = -0.0065f0 .+ 0.0003f0 .* annual_cycle
        vfal[:] = clamp.(0.58f0 .+ 0.18f0 .* annual_cycle, 0.12f0, 0.95f0)
        vslhf[:] = Float32(-9.0e3) .+ Float32(1.5e3) .* annual_cycle
        vsshf[:] = Float32(7.0e3) .+ Float32(1.3e3) .* shoulder_cycle
        vssrd[:] = max.(0.0f0, Float32(1.3e4) .+ Float32(2.7e3) .* annual_cycle)
        vstr[:] = Float32(-6.2e3) .+ Float32(1.0e3) .* shoulder_cycle

        ds.attrib["climate_source"] = "ERA5 CDS"
        ds.attrib["climate_frequency"] = "monthly"
        ds.attrib["ref_hgt"] = Float32(2500.0)
    finally
        close(ds)
    end
end

function _ensure_era5_file_for_tests(rgi_path::String)
    monthly = joinpath(rgi_path, "climate_historical_monthly_ERA5.nc")
    if !isfile(monthly)
        _write_synthetic_era5_monthly(monthly)
    end
end

function climate_downscale(; save_refs::Bool = false, climate_data_source::Symbol = :W5E5)
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]

    source_suffix = climate_data_source == :W5E5 ? "" : "_$(String(climate_data_source))"
    climate_ref_dir = joinpath(Sleipnir.root_dir, "test/data/climate")
    climate_step_ref_path = joinpath(climate_ref_dir, "climate_step$(source_suffix).jld2")
    climate_step_period_ref_path = joinpath(
        climate_ref_dir, "climate_step_period$(source_suffix).jld2")
    climate_2D_step_ref_path = joinpath(
        climate_ref_dir, "climate_2D_step$(source_suffix).jld2")

    tmpdir = nothing
    try
        if climate_data_source == :ERA5
            # Write the synthetic ERA5 fixture into a dedicated temp directory so
            # that the real prepro directory is never modified by tests.
            tmpdir = mktempdir()
            real_glacier_dir = joinpath(Sleipnir.prepro_dir, rgi_paths[rgi_ids[1]])
            tmp_glacier_dir = joinpath(tmpdir, rgi_ids[1])
            mkpath(tmp_glacier_dir)
            for f in readdir(real_glacier_dir; join = true)
                fname = basename(f)
                if !startswith(fname, "climate_") && !startswith(fname, "raw_climate_")
                    cp(f, joinpath(tmp_glacier_dir, fname))
                end
            end
            _ensure_era5_file_for_tests(tmp_glacier_dir)
            rgi_paths = Dict{String, String}(rgi_ids[1] => tmp_glacier_dir)
        end

        params = Parameters(
            simulation = SimulationParameters(
            use_velocities = false,
            use_glathida_data = true,
            multiprocessing = false,
            workers = 1,
            working_dir = Sleipnir.root_dir,
            test_mode = true,
            rgi_paths = rgi_paths,
            climate_data_source = climate_data_source
        )
        )

        rgi_path = joinpath(Sleipnir.prepro_dir, params.simulation.rgi_paths[rgi_ids[1]])
        if climate_data_source == :W5E5
            @test isfile(joinpath(rgi_path, "climate_historical_daily_W5E5.nc"))
        else
            @test isfile(joinpath(rgi_path, "climate_historical_monthly_ERA5.nc"))
        end

        glacier = initialize_glaciers(rgi_ids, params)[1]
        @test glacier.climate.climate_data_source == climate_data_source

        step = 1/12
        t = 2011.0
        period = partial_year(Day, t - step):Day(1):partial_year(Day, t)

        # Perform climate downscaling
        climate_step = get_cumulative_climate(glacier.climate.raw_climate)
        get_cumulative_climate!(glacier.climate, t, step)
        climate_2D_step = downscale_2D_climate(
            glacier.climate.climate_step, glacier.S, glacier.Coords)
        downscale_2D_climate!(glacier)

        JET.@test_opt broken=true target_modules=(Sleipnir,) get_cumulative_climate(glacier.climate.raw_climate)
        JET.@test_opt broken=false target_modules=(Sleipnir,) get_cumulative_climate!(
            glacier.climate, t, step)
        JET.@test_opt broken=false target_modules=(Sleipnir,) downscale_2D_climate!(glacier)

        if save_refs
            jldsave(climate_step_ref_path; climate_step)
            climate_step_period = glacier.climate.climate_step
            jldsave(climate_step_period_ref_path; climate_step_period)
            jldsave(climate_2D_step_ref_path; climate_2D_step)
        end

        @test isfile(climate_step_ref_path)
        @test isfile(climate_step_period_ref_path)
        @test isfile(climate_2D_step_ref_path)

        climate_step_ref = load(climate_step_ref_path)["climate_step"]
        climate_step_period_ref = load(climate_step_period_ref_path)["climate_step_period"]
        climate_2D_step_ref = load(climate_2D_step_ref_path)["climate_2D_step"]

        @test climate_step == climate_step_ref
        @test glacier.climate.climate_step == climate_step_period_ref
        @test climate_2D_step == climate_2D_step_ref
        @test glacier.climate.climate_2D_step == climate_2D_step_ref
    finally
        isnothing(tmpdir) || rm(tmpdir; recursive = true, force = true)
    end
end

function dummy_climate()
    climate = Sleipnir.DummyClimate2D(
        longterm_temps_scalar = [-2.0], longterm_temps_gridded = [-2.0 -1.5; -1.0 -0.5])
end
