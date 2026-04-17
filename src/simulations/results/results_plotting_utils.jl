###############################################################################
# Plotting utilities entry point
# Includes glacier, gridded-data, and video plotting sub-modules.
###############################################################################

export plot_glacier, plot_glacier_heatmaps, plot_glacier_quivers,
       plot_glacier_difference_evolution, plot_glacier_statistics_evolution,
       plot_glacier_integrated_volume, plot_glacier_dem
export plot_gridded_data, accumulate_gridded_data, plot_cumulative_gridded_data,
       plot_cumulative_mb, save_figure

using CairoMakie: Axis

include("results_plotting_glacier_utils.jl")
include("results_plotting_gridded_utils.jl")
include("results_plotting_video_utils.jl")
