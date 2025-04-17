export MeanDateVelocityMapping, IntegratedTrajectoryMapping

"""
    VelocityMapping

Abstract type representing the mapping to use in order to map the ice velocity
products onto the glacier grid.
It contains all needed information to build both the spatial projection, and how to
interpolate the data in time.
"""
abstract type VelocityMapping end

"""
    MeanDateVelocityMapping <: VelocityMapping

Mean date velocity mapping. It is the most simple mapping one can build and it
consists in taking the 2D vector field of ice velocity associated to a given mean
date and compare it to the instantaneous ice surface velocity obtained from the ice
flow model. It is valid only for ice surface velocities estimated from short time
windows since the velocity can vary within this time window.

# Fields
- `spatialInterp::String`: The spatial interpolation to use to map the ice surface
    velocity grid to the glacier grid. For the moment only "nearest" is supported.
"""
@kwdef struct MeanDateVelocityMapping <: VelocityMapping
    spatialInterp::String = "nearest"
end

"""
    IntegratedTrajectoryMapping <: VelocityMapping

Integrated trajectory mapping. This mapping is closer to reality as it consists in
integrating over time the instantaneous ice surface velocities along ice flow
trajectories in a Lagrangian way. This integrated velocity is then compared to the
velocity of the datacube. It has not been implemented yet but its computational cost
will likely be expensive.

# Fields
- `spatialInterp::String`: The spatial interpolation to use to map the ice surface
    velocity grid to the glacier grid. For the moment only "nearest" is supported.
"""
@kwdef struct IntegratedTrajectoryMapping <: VelocityMapping
    spatialInterp::String = "nearest"
end

"""
    grid(
        glacier::G,
        latitudes::Vector{F},
        longitudes::Vector{F},
        vx::Union{FileArray, Array{Union{Missing, F}, 3}},
        vy::Union{FileArray, Array{Union{Missing, F}, 3}},
        mapping::VM
    ) where {
        G <: AbstractGlacier,
        F <: AbstractFloat,
        VM <: VelocityMapping,
        FileArray <: Rasters.FileArray
    }

Grid velocity data onto the glacier grid following the prescribed mapping.
This function maps the 3 dimensional surface velocities to the glacier grid.
The provided surface velocities can be a `Rasters.FileArray` which happens when the
`RasterStack` is instantiated in lazy mode. In this situation, only the smallest
cube that contains all the needed data to construct the mapping is read from disk.
The returned velocity variables have shape `(nTimes, nx, ny)` where `nTimes` is the
number of time steps and `(nx, ny)` is the size of the glacier grid.

Arguments:
- `glacier::G`: Glacier instance which determines the glacier on which the
    velocities are projected onto.
- `latitudes::Vector{F}`: Vector of latitude values of the original surface
    velocity grid.
- `longitudes::Vector{F}`: Vector of longitude values of the original surface
    velocity grid.
- `vx::Union{FileArray, Array{Union{Missing, F}, 3}}`: X component of the original
    surface velocities. It can be either a `Rasters.FileArray` if the datacube is
    read in lazy mode, or a plain 3 dimensional array.
- `vy::Union{FileArray, Array{Union{Missing, F}, 3}}`: Y component of the original
    surface velocities. It can be either a `Rasters.FileArray` if the datacube is
    read in lazy mode, or a plain 3 dimensional array.
- `mapping::VM`: Mapping to use.

Returns:
- `xG`: A vector that gives the x coordinates of the glacier grid.
- `yG`: A vector that gives the y coordinates of the glacier grid.
- `vxG`: A 3 dimensional array of the x component of the velocity gridded onto the
    glacier grid.
- `vyG`: A 3 dimensional array of the y component of the velocity gridded onto the
    glacier grid.
"""
function grid(
    glacier::G,
    latitudes::Vector{F},
    longitudes::Vector{F},
    vx::Union{FileArray, Array{Union{Missing, F}, 3}},
    vy::Union{FileArray, Array{Union{Missing, F}, 3}},
    mapping::VM
) where {
    G <: AbstractGlacier,
    F <: AbstractFloat,
    VM <: VelocityMapping,
    FileArray <: Rasters.FileArray
}
    params_projection = glacier.params_projection
    transformReverse(lat,lon) = ReverseUTMercator(
        lat, lon;
        k=params_projection["k"],
        cenlon=params_projection["lon_0"], cenlat=params_projection["lat_0"],
        x0=params_projection["x_0"], y0=params_projection["y_0"]
    )
    xG = map(x -> x.x.val, transformReverse.(glacier.cenlat, glacier.Coords["lon"]))
    yG = map(x -> x.y.val, transformReverse.(glacier.Coords["lat"], glacier.cenlon))
    cenlatVel = (latitudes[end]+latitudes[begin])/2
    cenlonVel = (longitudes[end]+longitudes[begin])/2
    xV = map(x -> x.x.val, transformReverse.(cenlatVel, longitudes))
    yV = map(x -> x.y.val, transformReverse.(latitudes, cenlonVel))

    velType = eltype(vx)
    vxG = zeros(velType, glacier.nx, glacier.ny, size(vx,3))
    vyG = zeros(velType, glacier.nx, glacier.ny, size(vx,3))

    if mapping.spatialInterp=="nearest"
        # We express each of the coordinates of the velocity grid as (ΔxV*ix+bx, ΔyV*iy+by)
        # While this is an approximation since the coordinates do not truly live on a grid because of the projection, this error does not exceed one meter for most glaciers. This allows us to avoid having to compute an argmin.
        ΔxV = mean(diff(xV))
        ΔyV = mean(diff(yV))
        bx = xV[1]-ΔxV
        by = yV[1]-ΔyV
        ix = collect(1:length(longitudes))
        iy = collect(1:length(latitudes))
        indx = Int.(round.((xG .- bx)./ΔxV))
        indy = Int.(round.((yG .- by)./ΔyV))

        # Lazy arrays need to be read by block, hence we read the smallest block of data that contains all the points we need
        @assert size(vx,1)>=maximum(indx) "It looks like the datacube doesn't cover the whole glacier grid on the x-axis, please check that you use the right datacube for glacier $(glacier.rgi_id)."
        @assert size(vx,2)>=maximum(indy) "It looks like the datacube doesn't cover the whole glacier grid on the y-axis, please check that you use the right datacube for glacier $(glacier.rgi_id)."
        block_vx = vx[minimum(indx):maximum(indx),minimum(indy):maximum(indy),:]
        block_vy = vy[minimum(indx):maximum(indx),minimum(indy):maximum(indy),:]

        # Assign to each point of the glacier grid the closest point on the grid of surface velocities
        shiftx = -minimum(indx)+1
        shifty = -minimum(indy)+1
        for ix in range(1,length(glacier.Coords["lon"])), iy in range(1,length(glacier.Coords["lat"]))
            vxG[ix,iy,begin:end] .= block_vx[indx[ix]+shiftx,indy[iy]+shifty,:]
            vyG[ix,iy,begin:end] .= block_vy[indx[ix]+shiftx,indy[iy]+shifty,:]
        end
    else
        @error "$(mapping.spatialInterp) spatial interpolation method is not implemented"
    end
    return xG, yG, vxG, vyG
end
