###############################################################
########  COORDINATE PROJECTION UTILITIES  ####################
###############################################################

"""
    parse_proj(proj::String)

Parses the string containing the information of the projection to filter for important information
"+proj=tmerc +lat_0=0 +lon_0=6.985 +k=0.9996 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
"""
function parse_proj(proj::String)
    res = Dict()
    ℓ = split(proj, (' ', '+', '='))
    ℓ = ℓ[ℓ .!= ""]
    for (i, key) in enumerate(ℓ)
        if key ∈ ["lat_0", "lon_0", "k", "x_0", "y_0", "zone"]
            res[key] = parse(Float64, ℓ[i + 1])
        end
    end
    return res
end

"""
    UTMercator(x::F, y::F; k=0.9996, cenlon=0.0, cenlat=0.0, x0=0.0, y0=0.0, zone::Union{Nothing, Int}=nothing, hemisphere=nothing) where {F <: AbstractFloat}

Transverse Mercator Projection.
This function reprojects northing/easting coordinates into latitude/longitude.

# Keyword arguments

    - `k`: scale factor of the projection
    - `cenlon`: Central longitude used in the projection
    - `cenlat`: Central latitude used in the projection
    - `x0`: Shift in easting
    - `y0`: Shift in northing
    - `zone` : Zone of the projection
    - `hemisphere`: Either :north or :south
"""
function UTMercator(x::F, y::F; k = 0.9996, cenlon = 0.0, cenlat = 0.0,
        x0 = 0.0, y0 = 0.0, zone::Union{Nothing, Int} = nothing,
        hemisphere = nothing) where {F <: AbstractFloat}
    if !isnothing(zone)
        @assert !isnothing(hemisphere) "When zone is provided, hemisphere should also be defined. It can be either :north or :south"
        projection = CoordRefSystems.utm(hemisphere, zone; datum = WGS84Latest)(x, y)
    else
        # Convert to right units
        lonₒ = cenlon * 1.0°
        latₒ = cenlat * 1.0°
        xₒ = x0 * 1.0m
        yₒ = y0 * 1.0m
        # Define shift in new coordinate system
        S = CoordRefSystems.Shift(; lonₒ, xₒ, yₒ)
        # Define custom projection
        projection = TransverseMercator{k, latₒ, WGS84Latest, S}(x, y)
    end

    return convert(LatLon, projection)
end

"""
    ReverseUTMercator(x::F, y::F; k=0.9996, cenlon=0.0, cenlat=0.0, x0=0.0, y0=0.0, zone::Union{Nothing, Int}=nothing, hemisphere=nothing) where {F <: AbstractFloat}

Transverse Mercator Projection.
This function reprojects latitude/longitude into northing/easting coordinates.

# Keyword arguments

    - `k`: scale factor of the projection
    - `cenlon`: Central longitude used in the projection
    - `cenlat`: Central latitude used in the projection
    - `x0`: Shift in easting
    - `y0`: Shift in northing
    - `zone` : Zone of the projection
    - `hemisphere`: Either :north or :south
"""
function ReverseUTMercator(lat::F, lon::F; k = 0.9996, cenlon = 0.0, cenlat = 0.0,
        x0 = 0.0, y0 = 0.0, zone::Union{Nothing, Int} = nothing,
        hemisphere = nothing) where {F <: AbstractFloat}
    if !isnothing(zone)
        @assert !isnothing(hemisphere) "When zone is provided, hemisphere should also be defined. It can be either :north or :south"
        projection = CoordRefSystems.utm(hemisphere, zone; datum = WGS84Latest)
    else
        # Convert to right units
        lonₒ = cenlon * 1.0°
        latₒ = cenlat * 1.0°
        xₒ = x0 * 1.0m
        yₒ = y0 * 1.0m
        # Define shift in new coordinate system
        S = CoordRefSystems.Shift(; lonₒ, xₒ, yₒ)
        # Define custom projection
        projection = TransverseMercator{k, latₒ, WGS84Latest, S}
    end

    latlon = LatLon(lat, lon)
    return convert(projection, latlon)
end
