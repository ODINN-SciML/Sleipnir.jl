###############################################################################
# Glacier video utilities (CairoMakie)
###############################################################################

export plot_glacier_vid

"""
Generate a thickness-evolution animation. The output format (e.g. MP4 or GIF) is
inferred from the extension of `pathVideo`, since Makie's `record` selects the
encoder from the file extension.
"""
function make_thickness_video(
        results::Results,
        glacier::Glacier2D,
        tspan,
        step,
        pathVideo::String;
        colormap::Symbol = :viridis,
        colorrange::Union{Tuple, Nothing} = nothing,
        framerate::Int = 24,
        baseTitle::String = "Ice thickness"
)
    H_plot = results.H

    if length(glacier.Coords["lat"]) > 0
        lat = glacier.Coords["lat"]
        lon = glacier.Coords["lon"]
    else
        lat = glacier.Δy .* collect(1:1:glacier.ny)
        lon = glacier.Δx .* collect(1:1:glacier.nx)
    end
    x = results.x
    y = results.y

    H = [ifelse.(h .== 0.0, NaN, h) for h in H_plot]

    fig = CairoMakie.Figure(; size = (600, 600))

    # Create an axis
    ax = CairoMakie.Axis(fig[1, 1], aspect = DataAspect())
    xlims!(ax, lon[begin], lon[end])
    ylims!(ax, lat[begin], lat[end])
    ax.xlabel = "Longitude (°)"
    ax.ylabel = "Latitude (°)"

    # Number of frames
    nFrames = size(H, 1)

    if isnothing(colorrange)
        maxH = maximum(maximum.(H_plot))
        colorrange = (0.0, maxH)
    end

    hm = CairoMakie.heatmap!(
        ax,
        lon,
        lat,
        reverseForHeatmap(H[1], x, y),
        colorrange = colorrange,
        colormap = colormap,
        nan_color = :transparent
    )
    CairoMakie.Colorbar(fig[1, 2], hm, label = "Ice thickness (m)")

    years = tspan[1] .+ step * collect(1:size(H, 1))

    # Function to update the heatmap for each frame
    function _update_heatmap(frame_nb)
        hm[3] = reverseForHeatmap(H[frame_nb], x, y)
        year = Int(floor(years[frame_nb]))
        ax.title = "$baseTitle (t = $year)"
    end

    # Record the animation
    record(fig, pathVideo, 1:nFrames; framerate = framerate) do frame
        _update_heatmap(frame)
    end
end

"""
    plot_glacier_vid(
        plot_type::String,
        results::Results,
        glacier::Glacier2D,
        tspan,
        step,
        pathVideo::String;
        framerate::Int=24,
        baseTitle::String=""
    )

Generate various types of videos for glacier data. For now only the evolution of the glacier ice thickness is supported. More types of visualizations will be added in the future.

# Arguments

  - `plot_type`: Type of plot to generate. Options are:

      + "thickness": Heatmap of the glacier thickness.

  - `results`: A result object containing the simulation results including ice
    thickness over time.

  - `glacier`: A glacier instance.

  - `tspan`: The simulation time span.

  - `step`: Time step to use to retrieve the results and generate the video.

  - `pathVideo`: Path of the output animation. The format is inferred from the file
    extension — e.g. `.mp4` for a video or `.gif` for an animated GIF.

# Optional Keyword Arguments

  - `framerate`: The framerate to use for the video generation.
  - `baseTitle`: The prefix to use in the title of the frames. In each frame it is
    concatenated with the value of the year in the form " (t=XXXX)".
"""
function plot_glacier_vid(
        plot_type::String,
        results::Results,
        glacier::Glacier2D,
        tspan,
        step,
        pathVideo::String;
        framerate::Int = 24,
        baseTitle::String = "Ice thickness"
)
    if plot_type == "thickness"
        make_thickness_video(results, glacier, tspan, step, pathVideo;
            framerate = framerate, baseTitle = baseTitle)
    else
        error("Invalid plot_type: $plot_type")
    end
end
