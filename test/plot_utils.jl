function glaciers2D_plots()
    rgi_paths = get_rgi_paths()
    rgi_ids = ["RGI60-07.00042"]
    rgi_paths = Dict(k => rgi_paths[k] for k in rgi_ids)

    params = Parameters(
        simulation=SimulationParameters(
            use_velocities=false,
            use_glathida_data=false,
            working_dir=Sleipnir.root_dir,
            test_mode=true,
            rgi_paths=rgi_paths
        )
    )
    glaciers = initialize_glaciers(rgi_ids, params)

    S = glaciers[1].S
    ifm = SimpleIceflowModel{Sleipnir.Float}(S)

    results = Results(
        glaciers[1],
        ifm;
        H = [abs.(randn(size(glaciers[1].H₀)...)),abs.(randn(size(glaciers[1].H₀)...))],
        Vx = [abs.(randn(size(glaciers[1].H₀)...)),abs.(randn(size(glaciers[1].H₀)...))],
        Vy = [abs.(randn(size(glaciers[1].H₀)...)),abs.(randn(size(glaciers[1].H₀)...))],
        Vx_ref = [abs.(randn(size(glaciers[1].H₀)...)),abs.(randn(size(glaciers[1].H₀)...))],
        Vy_ref = [abs.(randn(size(glaciers[1].H₀)...)),abs.(randn(size(glaciers[1].H₀)...))],
    )

    # Test execution
    @testset "plot_glacier tests" begin
        @testset "Heatmaps" begin
            plot_glacier(results, "heatmaps", [:H,:B]; plotContour=true)
        end

        @testset "Quivers" begin
            plot = plot_glacier(results, "quivers", [:V_ref, :V])
        end

        @testset "Statistics Evolution" begin
            plot_glacier(results, "evolution statistics", [:H], tspan=(2010.0,2015.0), metrics=["average","std","max","median","min"])
        end

        @testset "Difference Evolution" begin
            plot_glacier(results, "evolution difference", [:H], tspan=(2010.0,2015.0), metrics=["difference","hist"], figsize=(800,600))
        end

        @testset "Integrated Volume" begin
            plot_glacier(results, "integrated volume", [:H], tspan=(2010.0,2015.0))
        end

        @testset "Bias" begin
            plot_glacier(results, "bias", [:B,:S])
        end
    end
end


# Fake iceflow model struct to test the plotting functions
struct IceflowModelTest{F <: AbstractFloat} <: AbstractModel
    S::Matrix{F}
end
function make_thickness_video_test()
    rgi_ids = ["RGI60-11.03646"]
    rgi_paths = get_rgi_paths()
    working_dir = joinpath(Sleipnir.root_dir, "test/data")

    δt = 1/12
    params = Parameters(
        simulation = SimulationParameters(
            use_MB = true,
            use_iceflow = true,
            use_velocities = true,
            use_glathida_data = false,
            tspan = (2014.0, 2015.0),
            step_MB = δt,
            working_dir = working_dir,
            multiprocessing = true,
            workers = 1,
            rgi_paths = rgi_paths,
            ice_thickness_source = "Farinotti19",
        ),
    )

    glaciers = initialize_glaciers(rgi_ids, params)

    nSteps = (params.simulation.tspan[2]-params.simulation.tspan[1])/δt
    timeSteps = params.simulation.tspan[1].+collect(0:nSteps).*δt
    H = [rand(glaciers[1].nx, glaciers[1].ny) for t in timeSteps]

    tempPath = mktempdir()*".mp4"

    ifm = IceflowModelTest{Float64}(glaciers[1].S)
    results = Results(glaciers[1], ifm; H=H)

    plot_glacier_vid("thickness", results, glaciers[1], params.simulation.tspan, δt, tempPath; baseTitle="Bossons glacier")
end
