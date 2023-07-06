
mutable struct Results{F <: AbstractFloat} 
    rgi_id::String
    H::Vector{Matrix{F}}
    S::Matrix{F}
    B::Matrix{F}
    V::Matrix{F}
    Vx::Matrix{F}
    Vy::Matrix{F}
end


function Results(glacier::AbstractGlacier, ifm::IF;
        rgi_id::String = glacier.rgi_id,
        H::Vector{Matrix{F}} = Matrix{F}([]),
        S::Matrix{F} = zeros(F, size(ifm.S)),
        B::Matrix{F} = zeros(F, size(ifm.B)),
        V::Matrix{F} = zeros(F, size(ifm.V)),
        Vx::Matrix{F} = zeros(F, size(ifm.Vx)),
        Vy::Matrix{F} = zeros(F, size(ifm.Vy))
            ) where {F <: AbstractFloat, IF <: AbstractModel}

    # Build the results struct based on input values
    results = Results(rgi_id, H, S, B,
                      V, Vx, Vy)

    return results
end
