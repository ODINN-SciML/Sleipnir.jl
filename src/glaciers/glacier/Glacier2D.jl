
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

`Glacier2D{F <: AbstractFloat, I <: Integer}`

# Fields
- `rgi_id::Union{String, Nothing}`: The RGI (Randolph Glacier Inventory) identifier for the glacier.
- `name::String`: The name of the glacier if available.
- `climate::Union{Climate2D, Nothing}`: The climate data associated with the glacier.
- `H₀::Union{Matrix{F}, Nothing}`: Initial ice thickness matrix.
- `H_glathida::Union{Matrix{F}, Nothing}`: Ice thickness matrix from the GLATHIDA dataset.
- `S::Union{Matrix{F}, Nothing}`: Surface elevation matrix.
- `B::Union{Matrix{F}, Nothing}`: Bedrock elevation matrix.
- `V::Union{Matrix{F}, Nothing}`: Ice velocity magnitude matrix.
- `Vx::Union{Matrix{F}, Nothing}`: Ice velocity in the x-direction matrix.
- `Vy::Union{Matrix{F}, Nothing}`: Ice velocity in the y-direction matrix.
- `A::Union{F, Nothing}`: Flow law parameter.
- `C::Union{F, Nothing}`: Sliding law parameter.
- `n::Union{F, Nothing}`: Flow law exponent.
- `slope::Union{Matrix{F}, Nothing}`: Surface slope matrix.
- `dist_border::Union{Matrix{F}, Nothing}`: Distance to the glacier border matrix.
- `Coords::Union{Dict{String, Vector{Float64}}, Nothing}`: Coordinates dictionary with keys as coordinate names and values as vectors of coordinates.
- `Δx::Union{F, Nothing}`: Grid spacing in the x-direction.
- `Δy::Union{F, Nothing}`: Grid spacing in the y-direction.
- `nx::Union{I, Nothing}`: Number of grid points in the x-direction.
- `ny::Union{I, Nothing}`: Number of grid points in the y-direction.
- `cenlon::Union{F, Nothing}`: Longitude of the glacier center.
- `cenlat::Union{F, Nothing}`: Latitude of the glacier center.
- `params_projection::Dict{String, Float64}`: Projection parameters that allows mapping the regional grid to global WGS84 coordinates.
"""
mutable struct Glacier2D{F <: AbstractFloat, I <: Integer} <: AbstractGlacier
    rgi_id::Union{String, Nothing}
    name::String
    climate::Union{Climate2D, Nothing}
    H₀::Union{Matrix{F}, Nothing}
    H_glathida::Union{Matrix{F}, Nothing}
    S::Union{Matrix{F}, Nothing}
    B::Union{Matrix{F}, Nothing}
    V::Union{Matrix{F}, Nothing}
    Vx::Union{Matrix{F}, Nothing}
    Vy::Union{Matrix{F}, Nothing}
    A::Union{F, Nothing}
    C::Union{F, Nothing}
    n::Union{F, Nothing}
    slope::Union{Matrix{F}, Nothing}
    dist_border::Union{Matrix{F}, Nothing}
    Coords::Union{Dict{String, Vector{Float64}}, Nothing}
    Δx::Union{F, Nothing}
    Δy::Union{F, Nothing}
    nx::Union{I, Nothing}
    ny::Union{I, Nothing}
    cenlon::Union{F, Nothing}
    cenlat::Union{F, Nothing}
    params_projection::Dict{String, Float64}
    # data::Union{Vector{DAT}, Nothing} # We ideally want this, but not clear how to specify concrete DAT type.
    data::Union{Vector, Nothing}
end

"""
Constructs a `Glacier2D` object with the given parameters, including default ones.

    Glacier2D(;
        rgi_id::Union{String, Nothing} = nothing,
        name::String = "",
        climate::Union{Climate2D, Nothing} = nothing,
        H₀::Union{Matrix{F}, Nothing} = nothing,
        H_glathida::Union{Matrix{F}, Nothing} = nothing,
        S::Union{Matrix{F}, Nothing} = nothing,
        B::Union{Matrix{F}, Nothing} = nothing,
        V::Union{Matrix{F}, Nothing} = nothing,
        Vx::Union{Matrix{F}, Nothing} = nothing,
        Vy::Union{Matrix{F}, Nothing} = nothing,
        A::Union{F, Nothing} = nothing,
        C::Union{F, Nothing} = nothing,
        n::Union{F, Nothing} = nothing,
        slope::Union{Matrix{F}, Nothing} = nothing,
        dist_border::Union{Matrix{F}, Nothing} = nothing,
        Coords::Union{Dict{String, Vector{Float64}}, Nothing} = nothing,
        Δx::Union{F, Nothing} = nothing,
        Δy::Union{F, Nothing} = nothing,
        nx::Union{I, Nothing} = nothing,
        ny::Union{I, Nothing} = nothing,
        cenlon::Union{F, Nothing} = nothing,
        cenlat::Union{F, Nothing} = nothing,
        params_projection::Dict{String, Float64} = Dict{String, Float64}(),
        data::Union{Vector, Nothing} = nothing,
    ) where {F <: AbstractFloat, I <: Integer}

# Arguments
- `rgi_id::Union{String, Nothing}`: The RGI identifier for the glacier.
- `name::String`: The name of the glacier if available.
- `climate::Union{Climate2D, Nothing}`: The climate data associated with the glacier.
- `H₀::Union{Matrix{F}, Nothing}`: Initial ice thickness matrix.
- `H_glathida::Union{Matrix{F}, Nothing}`: Ice thickness matrix from GLATHIDA.
- `S::Union{Matrix{F}, Nothing}`: Surface elevation matrix.
- `B::Union{Matrix{F}, Nothing}`: Bed elevation matrix.
- `V::Union{Matrix{F}, Nothing}`: Ice velocity magnitude matrix.
- `Vx::Union{Matrix{F}, Nothing}`: Ice velocity in the x-direction matrix.
- `Vy::Union{Matrix{F}, Nothing}`: Ice velocity in the y-direction matrix.
- `A::Union{F, Nothing}`: Flow law parameter.
- `C::Union{F, Nothing}`: Sliding law parameter.
- `n::Union{F, Nothing}`: Flow law exponent.
- `slope::Union{Matrix{F}, Nothing}`: Slope matrix.
- `dist_border::Union{Matrix{F}, Nothing}`: Distance to border matrix.
- `Coords::Union{Dict{String, Vector{Float64}}, Nothing}`: Coordinates dictionary with keys "lon" and "lat".
- `Δx::Union{F, Nothing}`: Grid spacing in the x-direction.
- `Δy::Union{F, Nothing}`: Grid spacing in the y-direction.
- `nx::Union{I, Nothing}`: Number of grid points in the x-direction.
- `ny::Union{I, Nothing}`: Number of grid points in the y-direction.
- `cenlon::Union{F, Nothing}`: Central longitude of the glacier.
- `cenlat::Union{F, Nothing}`: Central latitude of the glacier.
- `params_projection::Dict{String, Float64}`: Projection parameters that allows mapping the regional grid to global WGS84 coordinates.

# Returns
- A `Glacier2D` object with the specified parameters.
"""
function Glacier2D(;
    rgi_id::Union{String, Nothing} = nothing,
    name::String = "",
    climate::Union{Climate2D, Nothing} = nothing,
    H₀::Union{Matrix{F}, Nothing} = nothing,
    H_glathida::Union{Matrix{F}, Nothing} = nothing,
    S::Union{Matrix{F}, Nothing} = nothing,
    B::Union{Matrix{F}, Nothing} = nothing,
    V::Union{Matrix{F}, Nothing} = nothing,
    Vx::Union{Matrix{F}, Nothing} = nothing,
    Vy::Union{Matrix{F}, Nothing} = nothing,
    A::Union{F, Nothing} = nothing,
    C::Union{F, Nothing} = nothing,
    n::Union{F, Nothing} = nothing,
    slope::Union{Matrix{F}, Nothing} = nothing,
    dist_border::Union{Matrix{F}, Nothing} = nothing,
    Coords::Union{Dict{String, Vector{Float64}}, Nothing} = nothing,
    Δx::Union{F, Nothing} = nothing,
    Δy::Union{F, Nothing} = nothing,
    nx::Union{I, Nothing} = nothing,
    ny::Union{I, Nothing} = nothing,
    cenlon::Union{F, Nothing} = nothing,
    cenlat::Union{F, Nothing} = nothing,
    params_projection::Dict{String, Float64} = Dict{String, Float64}(),
    # data::Union{Vector{DAT}, Nothing} = nothing,
    data::Union{Vector, Nothing} = nothing,
) where {F <: AbstractFloat, I <: Integer}

    # Define default float and integer type for constructor
    ft = Sleipnir.Float
    it = Sleipnir.Int

    return Glacier2D{ft,it}(rgi_id, name, climate, H₀, H_glathida, S, B, V, Vx, Vy, A, C, n, slope, dist_border, Coords, Δx, Δy, nx, ny, cenlon, cenlat, params_projection, data)
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
                                      a.params_projection == b.params_projection


Base.:(≈)(a::Glacier2D, b::Glacier2D) = a.rgi_id == b.rgi_id && a.name == b.name &&
                                        a.climate == b.climate &&
                                        safe_approx(a.H₀, b.H₀) && safe_approx(a.H_glathida, b.H_glathida) &&
                                        safe_approx(a.S, b.S) && safe_approx(a.B, b.B) && safe_approx(a.V, b.V) &&
                                        safe_approx(a.A, b.A) && safe_approx(a.C, b.C) && safe_approx(a.n, b.n) &&
                                        isapprox(a.slope, b.slope; rtol=1e-3) && safe_approx(a.dist_border, b.dist_border) &&
                                        safe_approx(a.Coords, b.Coords) && safe_approx(a.Δx, b.Δx) && safe_approx(a.Δy, b.Δy) &&
                                        safe_approx(a.nx, b.nx) && safe_approx(a.ny, b.ny) &&
                                        safe_approx(a.cenlon, b.cenlon) && safe_approx(a.cenlat, b.cenlat) &&
                                        safe_approx(a.params_projection, b.params_projection)

# Display setup
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
    if isnothing(glacier.H_glathida)
        printstyled("   w/o";color=:red)
    else
        printstyled("   w/";color=:blue)
    end
    println(" glathida elevation")

    if !isnothing(glacier.S) & !isnothing(glacier.H₀)
        print("Min,mean,max bedrock elevation S : ")
        printstyled("$(round(minimum(glacier.S[glacier.H₀.>0]);digits=1)) $(round(mean(glacier.S[glacier.H₀.>0]);digits=1)) $(round(maximum(glacier.S[glacier.H₀.>0]);digits=1))\n";color=:blue)
        print("Mean,max ice thickness H₀ : ")
        printstyled("$(round(mean(glacier.H₀[glacier.H₀.>0]);digits=1)) $(round(maximum(glacier.H₀[glacier.H₀.>0]);digits=1))\n";color=:blue)
    end

    print("A= ")
    printstyled(@sprintf("%.3e", glacier.A); color=:blue)
    print("  C= ")
    printstyled(glacier.C; color=:blue)
    print("  n= ")
    printstyled(glacier.n; color=:blue)
end

# Vectorial form
# If you don't understand what's going on below, refer to https://discourse.julialang.org/t/improving-doc-for-display-print-show-repr/69124/3
function Base.show(io::IO, ::MIME"text/plain", glaciers::Vector{G}) where {G <: AbstractGlacier}
    return Base.show(io, glaciers)
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


include("glacier2D_utils.jl")
include("../climate/climate2D_utils.jl")