
export Glacier2D, Climate2D

"""
    AbstractGlacier

An abstract type representing a glacier. This serves as a base type for different glacier implementations in the `Sleipnir` package.
"""
abstract type AbstractGlacier end

include("../climate/Climate2D.jl")

"""
A mutable struct representing a 2D glacier. Notice that all fields can be empty by
providing `nothing` as the default value.

/!\\ WARNING /!\\ `Glacier` objects should not be constructed
manually, but rather through the `initialize_glaciers` function.

`Glacier2D{F <: AbstractFloat, I <: Integer, CLIM <: Climate2D, THICKDATA <: Union{ThicknessData, Nothing}, SURFVELDATA <: Union{SurfaceVelocityData, Nothing}}`

# Fields
- `rgi_id::String`: The RGI (Randolph Glacier Inventory) identifier for the glacier.
- `name::String`: The name of the glacier if available.
- `climate::CLIM`: The climate data associated with the glacier.
- `H₀::Matrix{F}`: Initial ice thickness matrix.
- `H_glathida::Matrix{F}`: Ice thickness matrix from the GLATHIDA dataset.
- `S::Matrix{F}`: Surface elevation matrix.
- `B::Matrix{F}`: Bedrock elevation matrix.
- `V::Matrix{F}`: Ice velocity magnitude matrix.
- `Vx::Matrix{F}`: Ice velocity in the x-direction matrix.
- `Vy::Matrix{F}`: Ice velocity in the y-direction matrix.
- `A::F`: Flow law parameter.
- `C::F`: Sliding law parameter.
- `n::F`: Flow law exponent.
- `slope::Matrix{F}`: Surface slope matrix.
- `dist_border::Matrix{F}`: Distance to the glacier border matrix.
- `Coords::Dict{String, Vector{Float64}}`: Coordinates dictionary with keys as coordinate names and values as vectors of coordinates.
- `Δx::F`: Grid spacing in the x-direction.
- `Δy::F`: Grid spacing in the y-direction.
- `nx::I`: Number of grid points in the x-direction.
- `ny::I`: Number of grid points in the y-direction.
- `cenlon::F`: Longitude of the glacier center.
- `cenlat::F`: Latitude of the glacier center.
- `params_projection::Dict{String, Float64}`: Projection parameters that allows mapping the regional grid to global WGS84 coordinates.
- `thicknessData::THICKDATA`: Thickness data structure that is used to store the reference values.
- `velocityData::SURFVELDATA`: Surface velocity data structure that is used to store the reference values.
"""
mutable struct Glacier2D{F <: AbstractFloat, I <: Integer, CLIM <: Climate2D, THICKDATA <: Union{ThicknessData, Nothing}, SURFVELDATA <: Union{SurfaceVelocityData, Nothing}} <: AbstractGlacier
    rgi_id::String
    name::String
    climate::CLIM
    H₀::Matrix{F}
    H_glathida::Matrix{F}
    S::Matrix{F}
    B::Matrix{F}
    V::Matrix{F}
    Vx::Matrix{F}
    Vy::Matrix{F}
    A::F
    C::F
    n::F
    slope::Matrix{F}
    dist_border::Matrix{F}
    Coords::Dict{String, Vector{Float64}}
    Δx::F
    Δy::F
    nx::I
    ny::I
    cenlon::F
    cenlat::F
    params_projection::Dict{String, Float64}
    thicknessData::THICKDATA
    velocityData::SURFVELDATA
end

"""
Constructs a `Glacier2D` object with the given parameters, including default ones.

    Glacier2D(;
        rgi_id::String = "",
        name::String = "",
        climate::Climate2D = nothing,
        H₀::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        H_glathida::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        S::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        B::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        V::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        Vx::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        Vy::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        A::F = 0.0,
        C::F = 0.0,
        n::F = 0.0,
        slope::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        dist_border::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
        Coords::Dict{String, Vector{Float64}} = Dict{String, Vector{Float64}}("lon" => [], "lat" => []),
        Δx::F = 0.0,
        Δy::F = 0.0,
        nx::I = 0,
        ny::I = 0,
        cenlon::F = NaN,
        cenlat::F = NaN,
        params_projection::Dict{String, Float64} = Dict{String, Float64}(),
        thicknessData::THICKDATA = nothing,
        velocityData::SURFVELDATA = nothing,
    ) where {
        F <: AbstractFloat,
        I <: Integer,
        THICKDATA <: Union{ThicknessData, Nothing},
        SURFVELDATA <: Union{SurfaceVelocityData, Nothing},
    }

# Arguments
- `rgi_id::String`: The RGI identifier for the glacier.
- `name::String`: The name of the glacier if available.
- `climate::Climate2D`: The climate data associated with the glacier.
- `H₀::Matrix{F}`: Initial ice thickness matrix.
- `H_glathida::Matrix{F}`: Ice thickness matrix from GLATHIDA.
- `S::Matrix{F}`: Surface elevation matrix.
- `B::Matrix{F}`: Bed elevation matrix.
- `V::Matrix{F}`: Ice velocity magnitude matrix.
- `Vx::Matrix{F}`: Ice velocity in the x-direction matrix.
- `Vy::Matrix{F}`: Ice velocity in the y-direction matrix.
- `A::F`: Flow law parameter.
- `C::F`: Sliding law parameter.
- `n::F`: Flow law exponent.
- `slope::Matrix{F}`: Slope matrix.
- `dist_border::Matrix{F}`: Distance to border matrix.
- `Coords::Dict{String, Vector{Float64}}`: Coordinates dictionary with keys "lon" and "lat".
- `Δx::F`: Grid spacing in the x-direction.
- `Δy::F`: Grid spacing in the y-direction.
- `nx::I`: Number of grid points in the x-direction.
- `ny::I`: Number of grid points in the y-direction.
- `cenlon::F`: Central longitude of the glacier.
- `cenlat::F`: Central latitude of the glacier.
- `params_projection::Dict{String, Float64}`: Projection parameters that allows mapping the regional grid to global WGS84 coordinates.
- `thicknessData::THICKDATA`: Thickness data structure that is used to store the reference values.
- `velocityData::SURFVELDATA`: Surface velocity data structure that is used to store the reference values.

# Returns
- A `Glacier2D` object with the specified parameters.
"""
function Glacier2D(;
    rgi_id::String = "",
    name::String = "",
    climate::Climate2D = nothing,
    H₀::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    H_glathida::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    S::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    B::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    V::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    Vx::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    Vy::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    A::F = 0.0,
    C::F = 0.0,
    n::F = 0.0,
    slope::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    dist_border::Matrix{F} = Matrix{Sleipnir.Float}([;;]),
    Coords::Dict{String, Vector{Float64}} = Dict{String, Vector{Float64}}("lon" => [], "lat" => []),
    Δx::F = 0.0,
    Δy::F = 0.0,
    nx::I = 0,
    ny::I = 0,
    cenlon::F = NaN,
    cenlat::F = NaN,
    params_projection::Dict{String, Float64} = Dict{String, Float64}(),
    thicknessData::THICKDATA = nothing,
    velocityData::SURFVELDATA = nothing,
) where {
    F <: AbstractFloat,
    I <: Integer,
    THICKDATA <: Union{ThicknessData, Nothing},
    SURFVELDATA <: Union{SurfaceVelocityData, Nothing},
}

    # Define default float and integer type for constructor
    ft = Sleipnir.Float
    it = Sleipnir.Int

    return Glacier2D{ft,it,typeof(climate),typeof(thicknessData),typeof(velocityData)}(
        rgi_id, name, climate, H₀, H_glathida,
        S, B, V, Vx, Vy, A, C, n,
        slope, dist_border, Coords,
        Δx, Δy, nx, ny,
        cenlon, cenlat, params_projection,
        thicknessData, velocityData,
    )
end

"""
Copies a `Glacier2D` object and updates the thickness and/or surface velocity data.

    Glacier2D(
        glacier::Glacier2D;
        thicknessData::Union{ThicknessData, Nothing} = nothing,
        velocityData::Union{SurfaceVelocityData, Nothing} = nothing,
    )

# Arguments
- `glacier::Glacier2D`: The original glacier struct.
- `thicknessData::Union{ThicknessData, Nothing}`: Thickness data structure that is used to store the reference values. Default is `nothing` which keeps the existing thickness data.
- `velocityData::Union{SurfaceVelocityData, Nothing}`: Surface velocity data structure that is used to store the reference values. Default is `nothing` which keeps the existing surface velocity data.

# Returns
- A `Glacier2D` object that is a copy of the original one with the thickness and/or surface velocity data updated.
"""
function Glacier2D(
    glacier::Glacier2D;
    thicknessData::Union{ThicknessData, Nothing} = nothing,
    velocityData::Union{SurfaceVelocityData, Nothing} = nothing,
)
    # Define default float and integer type for constructor
    ft = Sleipnir.Float
    it = Sleipnir.Int
    thicknessData = isnothing(thicknessData) ? glacier.thicknessData : thicknessData
    velocityData = isnothing(velocityData) ? glacier.velocityData : velocityData
    return Glacier2D{ft, it,typeof(glacier.climate),typeof(thicknessData),typeof(velocityData)}(
        glacier.rgi_id, glacier.name, glacier.climate, glacier.H₀, glacier.H_glathida,
        glacier.S, glacier.B, glacier.V, glacier.Vx, glacier.Vy, glacier.A, glacier.C, glacier.n,
        glacier.slope, glacier.dist_border, glacier.Coords,
        glacier.Δx, glacier.Δy, glacier.nx, glacier.ny,
        glacier.cenlon, glacier.cenlat, glacier.params_projection,
        thicknessData, velocityData,
    )
end

###############################################
################### UTILS #####################
###############################################


Base.:(==)(a::Glacier2D, b::Glacier2D) = a.rgi_id == b.rgi_id && a.name == b.name &&
                                      a.climate == b.climate &&
                                      a.H₀ == b.H₀ && a.H_glathida == b.H_glathida && a.S == b.S && a.B == b.B && a.V == b.V &&
                                      a.A == b.A && a.C == b.C && a.n == b.n &&
                                      a.slope == b.slope && a.dist_border == b.dist_border &&
                                      a.Coords == b.Coords && a.Δx == b.Δx && a.Δy == b.Δy && a.nx == b.nx && a.ny == b.ny &&
                                      a.cenlon == b.cenlon && a.cenlat == b.cenlat &&
                                      a.params_projection == b.params_projection &&
                                      a.thicknessData == b.thicknessData && a.velocityData == b.velocityData


Base.:(≈)(a::Glacier2D, b::Glacier2D) = a.rgi_id == b.rgi_id && a.name == b.name &&
                                        a.climate == b.climate &&
                                        safe_approx(a.H₀, b.H₀) && safe_approx(a.H_glathida, b.H_glathida) &&
                                        safe_approx(a.S, b.S) && safe_approx(a.B, b.B) && safe_approx(a.V, b.V) &&
                                        a.A == b.A && a.C == b.C && a.n == b.n &&
                                        isapprox(a.slope, b.slope; rtol=1e-3) && safe_approx(a.dist_border, b.dist_border) &&
                                        safe_approx(a.Coords, b.Coords) && a.Δx == b.Δx && a.Δy == b.Δy &&
                                        a.nx == b.nx && a.ny == b.ny &&
                                        safe_approx(a.cenlon, b.cenlon) && safe_approx(a.cenlat, b.cenlat) &&
                                        safe_approx(a.params_projection, b.params_projection) &&
                                        safe_approx(a.thicknessData, b.thicknessData) && safe_approx(a.velocityData, b.velocityData)

diffToDict(a::Glacier2D, b::Glacier2D) = Dict{Symbol, Bool}(
    :rgi_id => a.rgi_id == b.rgi_id,
    :name => a.name == b.name,
    :climate => a.climate == b.climate,
    :H₀ => a.H₀ == b.H₀,
    :H_glathida => a.H_glathida == b.H_glathida,
    :S => a.S == b.S,
    :B => a.B == b.B,
    :V => a.V == b.V,
    :A => a.A == b.A,
    :C => a.C == b.C,
    :n => a.n == b.n,
    :slope => a.slope == b.slope,
    :dist_border => a.dist_border == b.dist_border,
    :Coords => a.Coords == b.Coords,
    :Δx => a.Δx == b.Δx,
    :nx => a.nx == b.nx,
    :ny => a.ny == b.ny,
    :cenlon => a.cenlon == b.cenlon,
    :cenlat => a.cenlat == b.cenlat,
    :params_projection => a.params_projection == b.params_projection,
    :thicknessData => a.thicknessData == b.thicknessData,
    :velocityData => a.velocityData == b.velocityData,
)

# Display setup
function Base.show(io::IO, type::MIME"text/plain", glacier::Glacier2D)
    Base.show(io, glacier)
end
function Base.show(io::IO, glacier::Glacier2D)
    if !isnothing(glacier.H₀)
        H=round.(255*glacier.H₀/maximum(glacier.H₀))
        display(Gray.(Int.(H'[end:-1:1,1:end])/255))
    end

    print("Glacier2D ")
    strName = glacier.name=="" ? "" : " ($(glacier.name))"
    printstyled("$(glacier.rgi_id)$(strName)";color=:yellow)
    print(" with a ")
    printstyled("$(glacier.nx)x$(glacier.ny)";color=:red)
    print(" grid ")
    printstyled("(Δx=$(glacier.Δx), Δy=$(glacier.Δy))";color=:red)
    println("")
    if !isnothing(glacier.cenlat) & !isnothing(glacier.cenlon)
        print("at position ")
        printstyled("($(round(glacier.cenlat;digits=6))°, $(round(glacier.cenlon;digits=6))°)";color=:green)
    else
        print("at undefined location")
    end
    if size(glacier.H_glathida) == (0, 0)
        printstyled("   w/o";color=:red)
    else
        printstyled("   w/";color=:blue)
    end
    println(" glathida elevation")

    print("Min,mean,max bedrock elevation S : ")
    printstyled("$(round(minimum(glacier.S[glacier.H₀.>0]);digits=1)) $(round(mean(glacier.S[glacier.H₀.>0]);digits=1)) $(round(maximum(glacier.S[glacier.H₀.>0]);digits=1))\n";color=:blue)
    print("Mean,max ice thickness H₀ : ")
    printstyled("$(round(mean(glacier.H₀[glacier.H₀.>0]);digits=1)) $(round(maximum(glacier.H₀[glacier.H₀.>0]);digits=1))\n";color=:blue)

    print("A= ")
    printstyled(@sprintf("%.3e", glacier.A); color=:blue)
    print("  C= ")
    printstyled(glacier.C; color=:blue)
    print("  n= ")
    printstyled(glacier.n; color=:blue)

    if isnothing(glacier.thicknessData)
        printstyled("\nw/o";color=:red)
    else
        printstyled("\nw/";color=:blue)
    end
    print(" thickness data  ")
    if isnothing(glacier.velocityData)
        printstyled("   w/o";color=:red)
    else
        printstyled("   w/";color=:blue)
    end
    print(" velocity data")
end

# Vectorial form
# If you don't understand what's going on below, refer to https://discourse.julialang.org/t/improving-doc-for-display-print-show-repr/69124/3
function Base.show(io::IO, type::MIME"text/plain", glaciers::Vector{G}) where {G <: AbstractGlacier}
    Base.show(io, glaciers)
end
function Base.show(io::IO, glaciers::Vector{G}) where {G <: AbstractGlacier}
    len = length(glaciers)
    print("$(len)-element $(typeof(glaciers))")
    try
        regions = counter([split(split(glacier.rgi_id, "-")[2], ".")[1] for glacier in glaciers])
        regionsFormatted = ["$(k[1]) (x$(k[2]))" for k in regions]
        println(" distributed over regions $(join(regionsFormatted, ", "))")
    catch
        println(" distributed over undefined regions")
    end
    if len>5
        print(join([glacier.rgi_id for glacier in glaciers[1:2]], " "))
        print(" ... ")
        println(join([glacier.rgi_id for glacier in glaciers[len-1:len]], " "))
    else
        println(join([glacier.rgi_id for glacier in glaciers], " "))
    end
end
