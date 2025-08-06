import Pkg
Pkg.activate(dirname(Base.current_project()))

using Sleipnir
using BenchmarkTools
using Logging
Logging.disable_logging(Logging.Info)
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10

println("# Performance benchmark")

rgi_paths = get_rgi_paths()
rgi_ids = ["RGI60-07.00042", "RGI60-07.00065"] # Use glaciers that have glathida data

params = Parameters(
    simulation=SimulationParameters(
        use_velocities=false,
        use_glathida_data=true,
        working_dir=Sleipnir.root_dir,
        test_mode=true,
        rgi_paths=rgi_paths
    )
)

println("## Benchmark of glaciers initialization")
t = @benchmark initialize_glaciers($rgi_ids, $params)
display(t)
println("")
