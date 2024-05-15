
mutable struct Results{F <: AbstractFloat} 
    rgi_id::String
    H::Vector{Matrix{F}}
    H_glathida::Union{Nothing, Vector{Matrix{F}}}
    S::Matrix{F}
    B::Matrix{F}
    V::Matrix{F}
    Vx::Matrix{F}
    Vy::Matrix{F}
    V_ref::Matrix{F}
    Vx_ref::Matrix{F}
    Vy_ref::Matrix{F}
    Δx::F             
    Δy::F  
    lon::F 
    lat::F
    θ::Union{Nothing, Vector{F}}
    loss::Union{Nothing, Vector{F}}          
end


function Results(glacier::G, ifm::IF;
        rgi_id::String = glacier.rgi_id,
        H::Vector{Matrix{F}} = Vector{Matrix{F}}([]),
        H_glathida::Union{Nothing, Vector{Matrix{F}}} = glacier.H_glathida,
        S::Matrix{F} = zeros(F, size(ifm.S)),
        B::Matrix{F} = zeros(F, size(ifm.B)),
        V::Matrix{F} = zeros(F, size(ifm.V)),
        Vx::Matrix{F} = zeros(F, size(ifm.Vx)),
        Vy::Matrix{F} = zeros(F, size(ifm.Vy)),
        V_ref::Matrix{F} = glacier.V,
        Vx_ref::Matrix{F} = glacier.Vx,
        Vy_ref::Matrix{F} = glacier.Vy,
        Δx::F = glacier.Δx,                
        Δy::F = glacier.Δy,   
        lon::F = glacier.cenlon,
        lat::F = glacier.cenlat,
        θ::Union{Nothing,Vector{F}} = nothing,
        loss::Union{Nothing,Vector{F}} = nothing                 
    ) where {G <: AbstractGlacier, F <: AbstractFloat, IF <: AbstractModel}

    # Build the results struct based on input values
    results = Results(rgi_id, H, H_glathida, S, B,
                      V, Vx, Vy, V_ref, Vx_ref, Vy_ref,
                      Δx, Δy,lon,lat,
                      θ, loss)                

    return results
end



