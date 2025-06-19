export MeanDateVelocityMapping, IntegratedTrajectoryMapping

"""
    VelocityMapping

Abstract type representing the mapping to use in order to map the ice velocity
products onto the glacier grid.
It contains all needed information to build both the spatial projection, and how to
interpolate the data in time.
"""
abstract type VelocityMapping end

"""
    MeanDateVelocityMapping <: VelocityMapping

Mean date velocity mapping. It is the most simple mapping one can build and it
consists in taking the 2D vector field of ice velocity associated to a given mean
date and compare it to the instantaneous ice surface velocity obtained from the ice
flow model. It is valid only for ice surface velocities estimated from short time
windows since the velocity can vary within this time window.

# Fields
- `spatialInterp::Symbol`: The spatial interpolation to use to map the ice surface
    velocity grid to the glacier grid. For the moment only `:nearest` is supported.
- `thresDate::Float`: Threshold in days to accept or reject a surface velocity
    snapshot. This allows selecting surface velocity data that is close enough to
    the current time step in the loss function.
"""
@kwdef struct MeanDateVelocityMapping <: VelocityMapping
    spatialInterp::Symbol = :nearest
    thresDate::Float = 60/365
end

"""
    IntegratedTrajectoryMapping <: VelocityMapping

Integrated trajectory mapping. This mapping is closer to reality as it consists in
integrating over time the instantaneous ice surface velocities along ice flow
trajectories in a Lagrangian way. This integrated velocity is then compared to the
velocity of the datacube. It has not been implemented yet but its computational cost
will likely be expensive.

# Fields
- `spatialInterp::Symbol`: The spatial interpolation to use to map the ice surface
    velocity grid to the glacier grid. For the moment only `:nearest` is supported.
"""
@kwdef struct IntegratedTrajectoryMapping <: VelocityMapping
    spatialInterp::Symbol = :nearest
end
