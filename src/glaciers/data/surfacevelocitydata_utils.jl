export initialize_surfacevelocitydata
export plot_timeseries, plot_count

"""

"""
function initialize_surfacevelocitydata(file::String; interp=false)

    # Date of first adquisition
    date1 = ncread(file, "date1")
    # Date of second adquisition
    date2 = ncread(file, "date2")

    # Middle date
    if !interp 
        date_mean = 0.5 .* date1 .+ 0.5 .* date2
        date_mean_offset = date1_offset = date2_offset = 0
    else 
        date_mean = ncread(file, "mid_date")
        date_mean_offset = datetime2julian(DateTime("2015-08-14")) - 2400000.5
        date1_offset = datetime2julian(DateTime("2015-07-30")) - 2400000.5
        date2_offset = datetime2julian(DateTime("2015-08-29")) - 2400000.5
    end

    # Convert dates from Modified Julian Days to datetipes
    date1 = mjd.(date1 .+ date1_offset)
    date2 = mjd.(date2 .+ date2_offset)
    date_mean = mjd.(date_mean .+ date_mean_offset)
    date_error = date2 .- date1

    # We convert from Date to DateTime
    date_mean = DateTime.(date_mean)
    date1 = DateTime.(date1)
    date2 = DateTime.(date2)

    # Read data from netcdf file
    x = ncread(file, "x")
    y = ncread(file, "y")

    # Velocity in the x direction (m/yr)
    vx = ncread(file, "vx")              
    vy = ncread(file, "vy")    

    # Compute absolute velocity
    vabs = (vx.^2 .+ vy.^2).^0.5          
    
    if !interp
        vx_error = ncread(file, "error_vx")
        vy_error = ncread(file, "error_vy")
        vabs_error = (vx_error.^2 .+ vy_error.^2 ).^0.5
    else         
        vx_error = nothing
        vy_error = nothing
        vabs_error = nothing
    end

    # Run some basic tests 
    nx, ny, ntimes = size(vx)
    @assert length(date1) == length(date2) == ntimes
    @assert nx == ny == 250

    # Spatial preprocessing
    proj_zone = ncgetatt(file, "mapping", "utm_zone_number")
    transform(X,Y) = Sleipnir.UTMercator(X, Y; zone=proj_zone, hemisphere=:north) 
    Coordinates = transform.(x, y)
    latitudes = map(x -> x.lat.val, Coordinates)
    longitudes = map(x -> x.lon.val, Coordinates)

    return SurfaceVelocityData(x=x, y=y, lat=latitudes, lon=longitudes, vx=vx, vy=vy, vabs=vabs, vx_error=vx_error, vy_error=vy_error, vabs_error=vabs_error, date=date_mean, date1=date1, date2=date2, date_error=date_error)
end


# Plot if ice surface velocity data 

"""

"""
function plot_timeseries(data::SurfaceVelocityData, idx::Int, idy::Int; ignore_zeros=false, saveas::Union{Nothing, String}=nothing)

    vabs = data.vabs[idx, idy, :]

    # Define begining of hydrologic year
    yrs = 2015:2021
    mealt_max = [DateTime("$(yr)-08-01") for yr in yrs]
    mealt_max_raw = Dates.datetime2julian.(mealt_max) .- 2400000.5
    mealt_max = Date.(mealt_max)

    date_raw = Dates.datetime2julian.(data.date) .- 2400000.5
    date1_raw = Dates.datetime2julian.(data.date1) .- 2400000.5
    date2_raw = Dates.datetime2julian.(data.date2) .- 2400000.5

    # Unfortunately, Plots does not support horizonal error bar in Date format, so 
    # we do this manually for now...

    if ignore_zeros
        vmax = 1.2 * maximum(vabs[vabs .> 0.0])
        vmin = 0.8 * minimum(vabs[vabs .> 0.0])
    else
        vmax = maximum(vabs)
        vmin = minimum(vabs)
    end

    _plot = Plots.scatter(date_raw, vabs, label="Velocity", yerr=data.vabs_error, ms=3, msw=0.5)

    for i in 1:length(date_raw)
        Plots.plot!(_plot, [date1_raw[i], date2_raw[i]], [vabs[i], vabs[i]], lw=0.5, lc=:black, legend=false)
    end 

    vline!(mealt_max_raw, label="August 1st")

    Plots.plot!(fontfamily="Computer Modern",
                title="Ice Surface Velocities",
                titlefontsize=18,
                tickfontsize=15,
                legendfontsize=15,
                guidefontsize=18,
                xlabel="Date",
                ylabel="Velocity (m/yr)",
                xticks=(mealt_max_raw, mealt_max),
                ylimits=(vmin,vmax),
                #xlimits=(10^(-4),10^(-1)),
                legend=false,
                margin= 10mm,
                size=(1400,500),
                dpi=600)

    if isnothing(saveas)
        return plot 
    else
        Plots.savefig(_plot, saveas)
    end
end

"""

"""
function plot_count(data::SurfaceVelocityData; saveas::Union{Nothing, String}=nothing)

    fig_count = Figure(size=(800, 600), axis=(;title="Counts"));

    v_count = mean((data.vabs .!== 0.0) .* (.!isnan.(data.vabs)), dims=3)[:,:,1]

    max_count = maximum(v_count)

    ax_ct, hm = heatmap(fig_count[1,1], v_count, colorrange=(0.0, max_count));
    Colorbar(fig_count[:, end+1], hm);     

    if isnothing(saveas)
        return fig_count
    else
        save(saveas, fig_count)  
    end
end