"""
    mutable struct Results{F <: AbstractFloat, I <: Int}

A mutable struct to store the results of simulations.

# Fields
- `rgi_id::String`: Identifier for the RGI (Randolph Glacier Inventory).
- `H::Vector{Matrix{F}}`: Vector of matrices representing glacier ice thickness `H` over time.
- `H_glathida::Union{Nothing, Vector{Matrix{F}}}`: Optional vector of matrices for Glathida ice thicknesses.
- `H_ref::Union{Nothing, Vector{Matrix{F}}}`: Reference data for ice thickness.
- `S::Matrix{F}`: Glacier surface altimetry.
- `B::Matrix{F}`: Glacier bedrock.
- `V::Matrix{F}`: Glacier ice surface velocities.
- `Vx::Matrix{F}`: x-component of the glacier ice surface velocity `V`.
- `Vy::Matrix{F}`: y-component of the glacier ice surface velocity `V`.
- `V_ref::Union{Nothing, Matrix{F}}`: Reference data for glacier ice surface velocities `V`.
- `Vx_ref::Union{Nothing, Matrix{F}}`: Reference data for the x-component of the glacier ice surface velocity `Vx`.
- `Vy_ref::Union{Nothing, Matrix{F}}`: Reference data for the y-component of the glacier ice surface velocity `Vy`.
- `Δx::F`: Grid spacing in the x-direction.
- `Δy::F`: Grid spacing in the y-direction.
- `lon::Union{Nothing, F}`: Optional longitude value.
- `lat::Union{Nothing, F}`: Optional latitude value.
- `nx::I`: Number of grid points in the x-direction.
- `ny::I`: Number of grid points in the y-direction.
- `tspan::Vector{F}`: Time span of the simulation.
- `θ::Union{Nothing, ComponentArray{F}}`: Machine learning model parameters.
- `loss::Union{Nothing, Vector{F}}` Vector with evolution of loss function.
"""
mutable struct Results{F <: AbstractFloat, I <: Int}
    rgi_id::String
    H::Vector{Matrix{F}}
    H_glathida::Union{Nothing, Vector{Matrix{F}}}
    H_ref::Union{Nothing, Vector{Matrix{F}}}
    S::Matrix{F}
    B::Matrix{F}
    V::Union{Nothing, Vector{Matrix{F}}}
    Vx::Union{Nothing, Vector{Matrix{F}}}
    Vy::Union{Nothing, Vector{Matrix{F}}}
    V_ref::Union{Nothing, Vector{Matrix{F}}}
    Vx_ref::Union{Nothing, Vector{Matrix{F}}}
    Vy_ref::Union{Nothing, Vector{Matrix{F}}}
    Δx::F
    Δy::F
    lon::Union{Nothing, F}
    lat::Union{Nothing, F}
    nx::I
    ny::I
    t::Union{Vector{F}, Nothing}
    tspan::Union{Tuple{F, F}, Nothing}
    θ::Union{Nothing, ComponentArray{F}}
    loss::Union{Nothing, Vector{F}}
end

Base.:(==)(a::Results, b::Results) = a.rgi_id == b.rgi_id && a.H == b.H &&
                                    a.H_glathida == b.H_glathida && a.H_ref == b.H_ref &&
                                    a.S == b.S && a.B == b.B &&
                                    a.V == b.V && a.Vx == b.Vx && a.Vy == b.Vy &&
                                    a.V_ref == b.V_ref && a.Vx_ref == b.Vx_ref && a.Vy_ref == b.Vy_ref &&
                                    a.Δx == b.Δx && a.Δy == b.Δy &&
                                    a.lon == b.lon && a.lat == b.lat &&
                                    a.nx == b.nx && a.ny == b.ny && a.t == b.t &&
                                    a.tspan == b.tspan && a.θ == b.θ && a.loss == b.loss

"""
    Results(glacier::G, ifm::IF;
        rgi_id::String = glacier.rgi_id,
        H::Union{Nothing, Vector{Matrix{F}}} = nothing,
        H_glathida::Union{Nothing, Vector{Matrix{F}}} = glacier.H_glathida,
        H_ref::Union{Nothing, Vector{Matrix{F}}} = nothing,
        S::Union{Nothing, Matrix{F}} = nothing,
        B::Union{Nothing, Matrix{F}} = nothing,
        V::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Vx::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Vy::Union{Nothing, Vector{Matrix{F}}} = nothing,
        V_ref::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Vx_ref::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Vy_ref::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Δx::F = glacier.Δx,
        Δy::F = glacier.Δy,
        lon::Union{Nothing, F} = glacier.cenlon,
        lat::Union{Nothing, F} = glacier.cenlat,
        nx::I = glacier.nx,
        ny::I = glacier.ny,
        t::Union{Vector{F}, Nothing} = nothing,
        tspan::Union{Tuple{F, F}, Nothing} = nothing,
        θ::Union{Nothing,ComponentArray{F}} = nothing,
        loss::Union{Nothing,Vector{F}} = nothing
    ) where {G <: AbstractGlacier, F <: AbstractFloat, IF <: AbstractModel, I <: Int}

Construct a `Results` object for a glacier simulation.

# Arguments
- `glacier::G`: The glacier object, subtype of `AbstractGlacier`.
- `ifm::IF`: The model object, subtype of `AbstractModel`.
- `rgi_id::String`: The RGI identifier for the glacier. Defaults to `glacier.rgi_id`.
- `H::Union{Nothing, Vector{Matrix{F}}}`: Ice thickness matrices. Defaults to nothing.
- `H_glathida::Union{Nothing, Vector{Matrix{F}}}`: Ice thickness from GlaThiDa. Defaults to `glacier.H_glathida`.
- `H_ref::Union{Nothing, Vector{Matrix{F}}}`: Reference ice thickness. Defaults to nothing.
- `S::Union{Nothing, Matrix{F}}`: Surface elevation matrix. Defaults to a zero matrix of the same size as `ifm.S`.
- `B::Union{Nothing, Matrix{F}}`: Bed elevation matrix. Defaults to a zero matrix of the same size as `glacier.B`.
- `V::Union{Nothing, Vector{Matrix{F}}}`: Velocity magnitude matrix. Defaults to nothing.
- `Vx::Union{Nothing, Vector{Matrix{F}}}`: Velocity in the x-direction matrix. Defaults to nothing.
- `Vy::Union{Nothing, Vector{Matrix{F}}}`: Velocity in the y-direction matrix. Defaults to nothing.
- `V_ref::Union{Nothing, Vector{Matrix{F}}}`: Reference velocity magnitude matrix. Defaults to nothing.
- `Vx_ref::Union{Nothing, Vector{Matrix{F}}}`: Reference velocity in the x-direction matrix. Defaults to nothing.
- `Vy_ref::Union{Nothing, Vector{Matrix{F}}}`: Reference velocity in the y-direction matrix. Defaults to nothing.
- `Δx::F`: Grid spacing in the x-direction. Defaults to `glacier.Δx`.
- `Δy::F`: Grid spacing in the y-direction. Defaults to `glacier.Δy`.
- `lon::Union{Nothing, F}`: Longitude of the glacier center. Defaults to `glacier.cenlon`.
- `lat::Union{Nothing, F}`: Latitude of the glacier center. Defaults to `glacier.cenlat`.
- `nx::I`: Number of grid points in the x-direction. Defaults to `glacier.nx`.
- `ny::I`: Number of grid points in the y-direction. Defaults to `glacier.ny`.
- `tspan::Tuple(F, F)`: Timespan of the simulation.
- `θ::Union{Nothing, ComponentArray{F}}`: Model parameters. Defaults to `nothing`.
- `loss::Union{Nothing, Vector{F}}`: Loss values. Defaults to `nothing`.

# Returns
- `results::Results`: A `Results` object containing the simulation results.
"""
function Results(glacier::G, ifm::IF;
        rgi_id::String = glacier.rgi_id,
        H::Union{Nothing, Vector{Matrix{F}}} = nothing,
        H_glathida::Union{Nothing, Vector{Matrix{F}}} = glacier.H_glathida,
        H_ref::Union{Nothing, Vector{Matrix{F}}} = nothing,
        S::Union{Nothing, Matrix{F}} = nothing,
        B::Union{Nothing, Matrix{F}} = nothing,
        V::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Vx::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Vy::Union{Nothing, Vector{Matrix{F}}} = nothing,
        V_ref::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Vx_ref::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Vy_ref::Union{Nothing, Vector{Matrix{F}}} = nothing,
        Δx::F = glacier.Δx,
        Δy::F = glacier.Δy,
        lon::Union{Nothing, F} = glacier.cenlon,
        lat::Union{Nothing, F} = glacier.cenlat,
        nx::I = glacier.nx,
        ny::I = glacier.ny,
        t::Union{Vector{F}, Nothing} = nothing,
        tspan::Union{Tuple{F, F}, Nothing} = nothing,
        θ::Union{Nothing,ComponentArray{F}} = nothing,
        loss::Union{Nothing,Vector{F}} = nothing
    ) where {G <: AbstractGlacier, F <: AbstractFloat, IF <: AbstractModel, I <: Int}

    if isnothing(H)
        H = Vector{Matrix{F}}([])
    end
    if isnothing(S)
        S = zeros(F, size(ifm.S))
    end
    if isnothing(B)
        B = zeros(F, size(glacier.B))
    end
    # Build the results struct based on input values
    results = Results{F, I}(rgi_id, H, H_glathida, H_ref, S, B,
                      V, Vx, Vy, V_ref, Vx_ref, Vy_ref,
                      Δx, Δy,lon,lat, nx, ny, t, tspan,
                      θ, loss)

    return results
end
