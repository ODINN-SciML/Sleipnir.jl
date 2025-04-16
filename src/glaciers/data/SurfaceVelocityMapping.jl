export MeanDateVelocityMapping, IntegratedVelocityMapping

abstract type VelocityMapping end

@kwdef struct MeanDateVelocityMapping <: VelocityMapping
    spatialInterp::String = "nearest"
end

@kwdef struct IntegratedVelocityMapping <: VelocityMapping
    spatialInterp::String = "nearest"
end

function grid(glacier::G, latitudes::Vector{F}, longitudes::Vector{F}, vx, vy, mapping::VM) where {G <: AbstractGlacier, F <: AbstractFloat, VM <: VelocityMapping}
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
