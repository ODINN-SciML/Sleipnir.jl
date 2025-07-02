export mapVelocity

"""
    mapVelocity(
        velocityMapping::MeanDateVelocityMapping,
        velocityData::SurfaceVelocityData,
        t::AbstractFloat,
    )

Retrieve the reference ice surface velocity for a given time step.
This mapping uses the nearest snap shot available within a time window whose length
is controlled by `velocityMapping.thresDate`.
If no snapshot is found in the time window of length 2*velocityMapping.thresDate,
the returned ice velocity components are empty matrices and the returned boolean flag `useVel` is set to false.

# Arguments:
- `velocityMapping::MeanDateVelocityMapping`: Mapping to map the reference ice velocity to a target time step `t`.
- `velocityData::SurfaceVelocityData`: Surface velocity data. This is usually an attribute of a glacier.
- `t::AbstractFloat`: Current time step.

# Returns
- `Vx_ref`: Matrix of the x-component of the ice surface velocity.
- `Vy_ref`: Matrix of the y-component of the ice surface velocity.
- `V_ref`: Matrix of the ice surface velocity magnitude.
- `useVel`: Boolean indicating whether the returned ice surface velocity can be used
    or not. The value of this boolean depends on the success of the ice surface
    velocity mapping at the current time step `t`.
"""
function mapVelocity(
    velocityMapping::MeanDateVelocityMapping,
    velocityData::SurfaceVelocityData,
    t::AbstractFloat,
)
    # TODO: precompute dates in float before simulation and store them in velocityData
    date_Vref = datetime_to_floatyear.(velocityData.date)
    # date1_Vref = datetime_to_floatyear.(velocityData.date1)
    # date2_Vref = datetime_to_floatyear.(velocityData.date2)
    useVel = false
    Vx_ref = Vy_ref = V_ref = Matrix{Sleipnir.Float}([;;])
    if length(date_Vref)>0
        errTime, indVel = findmin(abs.(date_Vref .- t))
        if errTime <= velocityMapping.thresDate
            # Use available data only if the mean date is close enough to the current time
            Vx_ref = velocityData.vx[indVel]
            Vy_ref = velocityData.vy[indVel]
            V_ref = velocityData.vabs[indVel]
            useVel = true
        end
    end
    return Vx_ref, Vy_ref, V_ref, useVel
end
