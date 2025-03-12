"""
    mutable struct Results{F <: AbstractFloat, I <: Int}

A mutable struct to store the results of simulations.

# Fields
- `rgi_id::String`: Identifier for the RGI (Randolph Glacier Inventory).
- `H::Vector{Matrix{F}}`: Vector of matrices representing glacier ice thickness `H` over time.
- `H_glathida::Union{Nothing, Vector{Matrix{F}}}`: Optional vector of matrices for Glathida ice thicknesses.
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
    S::Matrix{F}
    B::Matrix{F}
    V::Matrix{F}
    Vx::Matrix{F}
    Vy::Matrix{F}
    V_ref::Union{Nothing, Matrix{F}}
    Vx_ref::Union{Nothing, Matrix{F}}
    Vy_ref::Union{Nothing, Matrix{F}}
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


"""
    Results(glacier::G, ifm::IF; rgi_id::String = glacier.rgi_id, H::Vector{Matrix{F}} = Vector{Matrix{F}}([]), 
            H_glathida::Union{Nothing, Vector{Matrix{F}}} = glacier.H_glathida, S::Matrix{F} = zeros(F, size(ifm.S)), 
            B::Matrix{F} = zeros(F, size(ifm.B)), V::Matrix{F} = zeros(F, size(ifm.V)), 
            Vx::Matrix{F} = zeros(F, size(ifm.Vx)), Vy::Matrix{F} = zeros(F, size(ifm.Vy)), 
            V_ref::Union{Nothing, Matrix{F}} = glacier.V, Vx_ref::Union{Nothing, Matrix{F}} = glacier.Vx, 
            Vy_ref::Union{Nothing, Matrix{F}} = glacier.Vy, Δx::F = glacier.Δx, Δy::F = glacier.Δy, 
            lon::Union{Nothing, F} = glacier.cenlon, lat::Union{Nothing, F} = glacier.cenlat, 
            nx::I = glacier.nx, ny::I = glacier.ny, θ::Union{Nothing, ComponentArray{F}} = nothing,
            loss::Union{Nothing, Vector{F}} = Nothing) where {G <: AbstractGlacier, F <: AbstractFloat, IF <: AbstractModel}

Construct a `Results` object for a glacier simulation.

# Arguments
- `glacier::G`: The glacier object, subtype of `AbstractGlacier`.
- `ifm::IF`: The model object, subtype of `AbstractModel`.
- `rgi_id::String`: The RGI identifier for the glacier. Defaults to `glacier.rgi_id`.
- `H::Vector{Matrix{F}}`: Ice thickness matrices. Defaults to an empty vector.
- `H_glathida::Union{Nothing, Vector{Matrix{F}}}`: Ice thickness from GlaThiDa. Defaults to `glacier.H_glathida`.
- `S::Matrix{F}`: Surface elevation matrix. Defaults to a zero matrix of the same size as `ifm.S`.
- `B::Matrix{F}`: Bed elevation matrix. Defaults to a zero matrix of the same size as `ifm.B`.
- `V::Matrix{F}`: Velocity magnitude matrix. Defaults to a zero matrix of the same size as `ifm.V`.
- `Vx::Matrix{F}`: Velocity in the x-direction matrix. Defaults to a zero matrix of the same size as `ifm.Vx`.
- `Vy::Matrix{F}`: Velocity in the y-direction matrix. Defaults to a zero matrix of the same size as `ifm.Vy`.
- `V_ref::Union{Nothing, Matrix{F}}`: Reference velocity magnitude matrix. Defaults to `glacier.V`.
- `Vx_ref::Union{Nothing, Matrix{F}}`: Reference velocity in the x-direction matrix. Defaults to `glacier.Vx`.
- `Vy_ref::Union{Nothing, Matrix{F}}`: Reference velocity in the y-direction matrix. Defaults to `glacier.Vy`.
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
        H::Vector{Matrix{F}} = Vector{Matrix{F}}([]),
        H_glathida::Union{Nothing, Vector{Matrix{F}}} = glacier.H_glathida,
        S::Matrix{F} = zeros(F, size(ifm.S)),
        B::Matrix{F} = zeros(F, size(ifm.B)),
        V::Matrix{F} = zeros(F, size(ifm.V)),
        Vx::Matrix{F} = zeros(F, size(ifm.Vx)),
        Vy::Matrix{F} = zeros(F, size(ifm.Vy)),
        V_ref::Union{Nothing, Matrix{F}} = glacier.V,
        Vx_ref::Union{Nothing, Matrix{F}} = glacier.Vx,
        Vy_ref::Union{Nothing, Matrix{F}} = glacier.Vy,
        Δx::F = glacier.Δx,
        Δy::F = glacier.Δy,
        lon::Union{Nothing, F}  = glacier.cenlon,
        lat::Union{Nothing, F}  = glacier.cenlat,
        nx::I = glacier.nx,
        ny::I = glacier.ny,
        t::Union{Vector{F}, Nothing} = nothing, 
        tspan::Union{Tuple{F, F}, Nothing} = nothing,
        θ::Union{Nothing,ComponentArray{F}} = nothing,
        loss::Union{Nothing,Vector{F}} = nothing                 
    ) where {G <: AbstractGlacier, F <: AbstractFloat, IF <: AbstractModel, I <: Int}

    # Build the results struct based on input values
    results = Results(rgi_id, H, H_glathida, S, B,
                      V, Vx, Vy, V_ref, Vx_ref, Vy_ref,
                      Δx, Δy,lon,lat, nx, ny, t, tspan,
                      θ, loss)                

    return results
end



