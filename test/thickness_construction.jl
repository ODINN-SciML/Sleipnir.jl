
function thickness_construction()
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]

    params = Parameters(
        simulation = SimulationParameters(
        use_velocities = false,
        use_glathida_data = true,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        rgi_paths = rgi_paths
    )
    )

    glacier = initialize_glaciers(rgi_ids, params)[1]
    t = [2010.0, 2010.1]
    H = [randn(glacier.nx, glacier.ny), randn(glacier.nx, glacier.ny)]

    thickness_data = ThicknessData(t, H)
    JET.@test_opt target_modules=(Sleipnir,) ThicknessData(t, H)

    glacier_with_thickness_data = Glacier2D(glacier, thicknessData = thickness_data)
    JET.@test_opt target_modules=(Sleipnir,) Glacier2D(
        glacier, thicknessData = thickness_data)

    @test thickness_data==glacier_with_thickness_data.thicknessData
    @test isnothing(thickness_data.σ)
    σ = [fill(10.0, size(Hk)) for Hk in H]
    @test ThicknessData(t = t, H = H, σ = σ) != thickness_data
    @test ThicknessData(t = t, H = H, σ = σ) ≈ ThicknessData(t = t, H = H, σ = deepcopy(σ))
end

"""
    glathida_sigma()

The uncertainty of a glathida cell is the radar error, 5 m + 3% of H, and 1.5 times the
spread of the measurements inside the cell. The cells without measurement get 0.
"""
function glathida_sigma()
    # Native grid of 4x4 cells: one cell with the measurements 90 and 110, one with 100
    # three times, one with a single measurement of 200
    grid = zeros(4, 4)
    count = zeros(4, 4)
    sumsq = zeros(4, 4)
    for (i, j, values) in ((1, 1, [90.0, 110.0, 100.0]), (2, 2, [100.0, 100.0, 100.0]),
        (4, 4, [200.0]))
        grid[i, j] = sum(values) / length(values)
        count[i, j] = length(values)
        sumsq[i, j] = sum(abs2, values)
    end
    σ = Sleipnir.glathida_cell_σ(grid, grid, count, sumsq, 1)
    radar(H) = 5.0 + 0.03H
    @test σ[1, 1] ≈ sqrt(radar(100.0)^2 + (1.5 * 10.0)^2)
    @test σ[2, 2] ≈ radar(100.0)
    # Fewer than 3 measurements: the median spread of the other cells (0 and 10)
    @test σ[4, 4] ≈ sqrt(radar(200.0)^2 + (1.5 * 5.0)^2)
    @test σ[3, 3] == 0.0

    # By blocks of 2: the measurements of the upper left block are 90, 110, 100 (x4)
    H = Sleipnir._coarsen_glathida(grid, 2)
    σ₂ = Sleipnir.glathida_cell_σ(H, grid, count, sumsq, 2)
    values = [90.0, 110.0, 100.0, 100.0, 100.0, 100.0]
    @test σ₂[1, 1] ≈ sqrt(radar(H[1, 1])^2 + (1.5 * std(values))^2)
    @test σ₂[1, 2] == 0.0

    # No cell with 3 measurements: the spread is not known, and it is 12 m and not 0, so that
    # the uncertainty is not the radar error alone
    grid₁ = zeros(4, 4)
    count₁ = zeros(4, 4)
    sumsq₁ = zeros(4, 4)
    grid₁[1, 1], count₁[1, 1], sumsq₁[1, 1] = 100.0, 1.0, 100.0^2
    σ₃ = Sleipnir.glathida_cell_σ(grid₁, grid₁, count₁, sumsq₁, 1)
    @test σ₃[1, 1] ≈ sqrt(radar(100.0)^2 + (1.5 * 12.0)^2)
end

function surface_construction()
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]

    params = Parameters(
        simulation = SimulationParameters(
        use_velocities = false,
        working_dir = Sleipnir.root_dir,
        test_mode = true,
        rgi_paths = rgi_paths
    )
    )

    glacier = initialize_glaciers(rgi_ids, params)[1]
    # A cell without value is NaN, and two equal maps with NaN are equal
    S₂ = glacier.S .+ 1.0
    S₂[1, 1] = NaN
    surface_data = SurfaceData(t = [2000.1, 2012.8], S = [glacier.S, S₂],
        σ = [zero(glacier.S), fill(2.0, size(S₂))])
    glacier_with_surface_data = Glacier2D(glacier, surfaceData = surface_data)

    @test glacier_with_surface_data.surfaceData == deepcopy(surface_data)
    @test glacier_with_surface_data.surfaceData ≈ deepcopy(surface_data)
    @test tdata(surface_data) == [2000.1, 2012.8]
    # Each DEM is an observation
    @test count(row -> row.kind == :S, observations(glacier_with_surface_data)) == 2
end
