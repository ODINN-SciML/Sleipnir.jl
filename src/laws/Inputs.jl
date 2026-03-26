export iTopoSlope, iTopoAspect

"""
    iTopoSlope{F<:AbstractFloat} <: AbstractInput

Input representing dynamic surface slope (degrees) computed from the current
glacier surface topography.
"""
struct iTopoSlope{F <: AbstractFloat} <: AbstractInput
    window_m::F
    function iTopoSlope(; window_m::F = Sleipnir.Float(200.0)) where {F <: AbstractFloat}
        new{F}(window_m)
    end
end

default_name(::iTopoSlope) = :slope

function get_input(inp::iTopoSlope, simulation, glacier_idx, t)
    _ = t
    glacier = simulation.glaciers[glacier_idx]
    return compute_surface_slope(glacier.S, glacier.Δx, glacier.Δy; window_m = inp.window_m)
end

"""
    iTopoAspect{F<:AbstractFloat} <: AbstractInput

Input representing dynamic surface aspect (degrees, [0, 360)) computed from the
current glacier surface topography.
"""
struct iTopoAspect{F <: AbstractFloat} <: AbstractInput
    window_m::F
    function iTopoAspect(; window_m::F = Sleipnir.Float(200.0)) where {F <: AbstractFloat}
        new{F}(window_m)
    end
end

default_name(::iTopoAspect) = :aspect

function get_input(inp::iTopoAspect, simulation, glacier_idx, t)
    _ = t
    glacier = simulation.glaciers[glacier_idx]
    return compute_surface_aspect(glacier.S, glacier.Δx, glacier.Δy; window_m = inp.window_m)
end
