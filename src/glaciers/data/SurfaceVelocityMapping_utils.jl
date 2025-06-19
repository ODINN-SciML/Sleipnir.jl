export mapVelocity

function mapVelocity(
    velocityMapping::MeanDateVelocityMapping,
    Vx_ref::Vector{Matrix{F}},
    Vy_ref::Vector{Matrix{F}},
    V_ref::Vector{Matrix{F}},
    date_Vref::Vector{F},
    t::F,
) where {F <: AbstractFloat}
    useVel = false
    Vxτ_ref = Vyτ_ref = Vτ_ref = Matrix{Sleipnir.Float}([;;])
    if length(date_Vref)>0
        errTime, indVel = findmin(abs.(date_Vref .- t))
        if errTime <= velocityMapping.thresDate
            # Use available data only if the mean date is close enough to the current time
            Vxτ_ref = Vx_ref[indVel]
            Vyτ_ref = Vy_ref[indVel]
            Vτ_ref = V_ref[indVel]
            useVel = true
        end
    end
    return Vxτ_ref, Vyτ_ref, Vτ_ref, useVel
end
