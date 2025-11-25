"""
    mutable struct Results{F <: AbstractFloat, I <: Integer}

A mutable struct to store the results of simulations.

# Fields

  - `rgi_id::String`: Identifier for the RGI (Randolph Glacier Inventory).
  - `H::Vector{Matrix{F}}`: Vector of matrices representing glacier ice thickness `H` over time.
  - `H_glathida::Matrix{F}`: Optional matrix for Glathida ice thicknesses.
  - `H_ref::Vector{Matrix{F}}`: Reference data for ice thickness.
  - `S::Matrix{F}`: Glacier surface altimetry.
  - `B::Matrix{F}`: Glacier bedrock.
  - `V::Matrix{F}`: Glacier ice surface velocities.
  - `Vx::Matrix{F}`: x-component of the glacier ice surface velocity `V`.
  - `Vy::Matrix{F}`: y-component of the glacier ice surface velocity `V`.
  - `V_ref::Matrix{F}`: Reference data for glacier ice surface velocities `V`.
  - `Vx_ref::Matrix{F}`: Reference data for the x-component of the glacier ice surface velocity `Vx`.
  - `Vy_ref::Matrix{F}`: Reference data for the y-component of the glacier ice surface velocity `Vy`.
  - `date_Vref::Vector{F}`: Date of velocity observation (mean of `date1` and `date2`).
  - `date1_Vref::Vector{F}`: First date of velocity acquisition.
  - `date2_Vref::Vector{F}`: Second date of velocity acquisition.
  - `Δx::F`: Grid spacing in the x-direction.
  - `Δy::F`: Grid spacing in the y-direction.
  - `lon::F`: Longitude of the glacier grid center.
  - `lat::F`: Latitude of the glacier grid center.
  - `nx::I`: Number of grid points in the x-direction.
  - `ny::I`: Number of grid points in the y-direction.
  - `tspan::Vector{F}`: Time span of the simulation.
"""
mutable struct Results{F <: AbstractFloat, I <: Integer}
    rgi_id::String
    H::Vector{Matrix{F}}
    H_glathida::Matrix{F}
    H_ref::Vector{Matrix{F}}
    S::Matrix{F}
    B::Matrix{F}
    x::Vector{F}
    y::Vector{F}
    V::Vector{Matrix{F}}
    Vx::Vector{Matrix{F}}
    Vy::Vector{Matrix{F}}
    V_ref::Vector{Matrix{F}}
    Vx_ref::Vector{Matrix{F}}
    Vy_ref::Vector{Matrix{F}}
    date_Vref::Vector{F}
    date1_Vref::Vector{F}
    date2_Vref::Vector{F}
    Δx::F
    Δy::F
    lon::F
    lat::F
    nx::I
    ny::I
    t::Vector{F}
    tspan::Tuple{F, F}
end

function Base.:(==)(a::Results, b::Results)
    (
        a.rgi_id == b.rgi_id && a.H == b.H &&
        a.H_glathida == b.H_glathida && a.H_ref == b.H_ref &&
        a.S == b.S && a.B == b.B && a.x == b.x && a.y == b.y &&
        a.V == b.V && a.Vx == b.Vx && a.Vy == b.Vy &&
        a.V_ref == b.V_ref && a.Vx_ref == b.Vx_ref && a.Vy_ref == b.Vy_ref &&
        a.date_Vref == b.date_Vref && a.date1_Vref == b.date1_Vref &&
        a.date2_Vref == b.date2_Vref &&
        a.Δx == b.Δx && a.Δy == b.Δy &&
        a.lon == b.lon && a.lat == b.lat &&
        a.nx == b.nx && a.ny == b.ny && a.t == b.t &&
        isequal(a.tspan, b.tspan)
    )
end

"""
    Results(glacier::G, ifm::IF;
        rgi_id::String = glacier.rgi_id,
        H::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        H_glathida::Matrix{F} = glacier.H_glathida,
        H_ref::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        S::Matrix{F} = zeros(Sleipnir.Float, size(ifm.S)),
        B::Matrix{F} = zeros(Sleipnir.Float, size(glacier.B)),
        V::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        Vx::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        Vy::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        V_ref::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        Vx_ref::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        Vy_ref::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        date_Vref::Vector{F} = Vector{Sleipnir.Float}([]),
        date1_Vref::Vector{F} = Vector{Sleipnir.Float}([]),
        date2_Vref::Vector{F} = Vector{Sleipnir.Float}([]),
        Δx::F = glacier.Δx,
        Δy::F = glacier.Δy,
        lon::F = glacier.cenlon,
        lat::F = glacier.cenlat,
        nx::I = glacier.nx,
        ny::I = glacier.ny,
        t::Vector{F} = Vector{Sleipnir.Float}([]),
        tspan::Tuple{F, F} = (NaN, NaN),
    ) where {G <: AbstractGlacier, F <: AbstractFloat, IF <: AbstractModel, I <: Integer}

Construct a `Results` object for a glacier simulation.

# Arguments

  - `glacier::G`: The glacier object, subtype of `AbstractGlacier`.
  - `ifm::IF`: The model object, subtype of `AbstractModel`.
  - `rgi_id::String`: The RGI identifier for the glacier. Defaults to `glacier.rgi_id`.
  - `H::Vector{Matrix{F}}`: Ice thickness matrices. Defaults to an empty vector.
  - `H_glathida::Matrix{F}`: Ice thickness from GlaThiDa. Defaults to `glacier.H_glathida`.
  - `H_ref::Vector{Matrix{F}}`: Reference ice thickness. Defaults to an empty vector.
  - `S::Matrix{F}`: Surface elevation matrix. Defaults to a zero matrix of the same size as `ifm.S`.
  - `B::Matrix{F}`: Bed elevation matrix. Defaults to a zero matrix of the same size as `glacier.B`.
  - `V::Vector{Matrix{F}}`: Velocity magnitude matrix. Defaults to an empty vector.
  - `Vx::Vector{Matrix{F}}`: Velocity in the x-direction matrix. Defaults to an empty vector.
  - `Vy::Vector{Matrix{F}}`: Velocity in the y-direction matrix. Defaults to an empty vector.
  - `V_ref::Vector{Matrix{F}}`: Reference velocity magnitude matrix. Defaults to an empty vector.
  - `Vx_ref::Vector{Matrix{F}}`: Reference velocity in the x-direction matrix. Defaults to an empty vector.
  - `Vy_ref::Vector{Matrix{F}}`: Reference velocity in the y-direction matrix. Defaults to an empty vector.
  - `date_Vref::Vector{F}`: Date of velocity observation (mean of `date1` and `date2`). Defaults to an empty vector.
  - `date1_Vref::Vector{F}`: First date of velocity acquisition. Defaults to an empty vector.
  - `date2_Vref::Vector{F}`: Second date of velocity acquisition. Defaults to an empty vector.
  - `Δx::F`: Grid spacing in the x-direction. Defaults to `glacier.Δx`.
  - `Δy::F`: Grid spacing in the y-direction. Defaults to `glacier.Δy`.
  - `lon::F`: Longitude of the glacier grid center. Defaults to `glacier.cenlon`.
  - `lat::F`: Latitude of the glacier grid center. Defaults to `glacier.cenlat`.
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
        H::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        H_glathida::Matrix{F} = glacier.H_glathida,
        H_ref::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        S::Matrix{F} = zeros(Sleipnir.Float, size(ifm.S)),
        B::Matrix{F} = zeros(Sleipnir.Float, size(glacier.B)),
        V::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        Vx::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        Vy::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        V_ref::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        Vx_ref::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        Vy_ref::Vector{Matrix{F}} = Vector{Matrix{Sleipnir.Float}}([[;;]]),
        date_Vref::Vector{F} = Vector{Sleipnir.Float}([]),
        date1_Vref::Vector{F} = Vector{Sleipnir.Float}([]),
        date2_Vref::Vector{F} = Vector{Sleipnir.Float}([]),
        Δx::F = glacier.Δx,
        Δy::F = glacier.Δy,
        lon::F = glacier.cenlon,
        lat::F = glacier.cenlat,
        nx::I = glacier.nx,
        ny::I = glacier.ny,
        t::Vector{F} = Vector{Sleipnir.Float}([]),
        tspan::Tuple{F, F} = (NaN, NaN)
) where {G <: AbstractGlacier, F <: AbstractFloat, IF <: AbstractModel, I <: Integer}
    x = glacier.Coords["lon"]
    y = glacier.Coords["lat"]

    # Build the results struct based on input values
    results = Results{Sleipnir.Float, Sleipnir.Int}(
        rgi_id, H, H_glathida, H_ref, S, B,
        x, y,
        V, Vx, Vy, V_ref, Vx_ref, Vy_ref,
        date_Vref, date1_Vref, date2_Vref,
        Δx, Δy, lon, lat, nx, ny, t, tspan
    )

    return results
end
