
"""
    stop_condition_tstops(u, t, integrator, tstops)

Check if the current time `t` is in the list of stop times `tstops`.

# Arguments

  - `u`: The current state of the system (not used in this function).
  - `t::AbstractFloat`: The current time.
  - `integrator`: The integrator object (not used in this function).
  - `tstops::Vector{<: AbstractFloat}`: A collection of times at which the integration should stop.

# Returns

  - `Bool`: `true` if `t` is in `tstops`, otherwise `false`.
"""
function stop_condition_tstops(
        u, t::AbstractFloat, integrator, tstops::Vector{<: AbstractFloat})
    t in tstops
end
