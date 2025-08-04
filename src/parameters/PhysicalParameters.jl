
"""
A structure representing physical parameters used in simulations.

    PhysicalParameters{F <: AbstractFloat}

# Fields
- `ρ::F`: Density of ice.
- `g::F`: Gravitational acceleration.
- `ϵ::F`: Regularization used in the square root of norms for AD numerical stability.
- `η₀::F`: Initial viscosity.
- `maxA::F`: Maximum value for `A` (Glen's coefficient).
- `minA::F`: Minimum value for `A` (Glen's coefficient).
- `maxTlaw::F`: Maximum value of Temperature used in simulations on fake law.
- `minTlaw::F`: Minimum value of Temperature used in simulations on fake law.
- `noise_A_magnitude::F`: Magnitude of noise in A.
"""
struct PhysicalParameters{F <: AbstractFloat} <: AbstractParameters
    ρ::F
    g::F
    ϵ::F
    η₀::F
    maxA::F
    minA::F
    maxTlaw::F
    minTlaw::F
    noise_A_magnitude::F
end

"""
Initialize the physical parameters of a model.

    PhysicalParameters(;
        ρ::Float64 = 900.0,
        g::Float64 = 9.81,
        ϵ::Float64 = 1e-10,
        η₀::F = 1.0,
        maxA::Float64 = 8e-17,
        minA::Float64 = 8.5e-20,
        maxTlaw::Float64 = 1.0,
        minTlaw::Float64 = -25.0,
        noise_A_magnitude::Float64 = 5e-18
        )

Keyword arguments
=================
    - `ρ`: Density of ice.
    - `g`: Gravitational acceleration.
    - `ϵ`: Regularization used in the square root of norms for AD numerical stability.
    - `η₀`: Initial viscosity.
    - `maxA`: Maximum value for `A` (Glen's coefficient).
    - `minA`: Minimum value for `A` (Glen's coefficient).
    - `maxTlaw`: Maximum value of Temperature used in simulations on fake law.
    - `minTlaw`: Minimum value of Temperature used in simulations on fake law.
    - `noise_A_magnitude`: Magnitude of noise added to A
"""
function PhysicalParameters(;
            ρ::F = 900.0,
            g::F = 9.81,
            ϵ::F = 1e-10,
            η₀::F = 1.0,
            maxA::F = 8e-17,
            minA::F = 8.5e-20,
            maxTlaw::F = 1.0,
            minTlaw::F = -25.0,
            noise_A_magnitude::F = 5e-18
            ) where {F <: AbstractFloat}
    # Build PhysicalParameters based on values
    ft = typeof(g)
    physical_parameters = PhysicalParameters{ft}(ρ, g, ϵ, η₀,
                                            maxA, minA,
                                            maxTlaw, minTlaw,
                                            noise_A_magnitude)

    return physical_parameters
end

Base.:(==)(a::PhysicalParameters, b::PhysicalParameters) = a.ρ == b.ρ && a.g == b.g &&
                                      a.ϵ == b.ϵ && a.η₀ == b.η₀ &&
                                      a.maxA == b.maxA && a.minA == b.minA && a.maxTlaw == b.maxTlaw &&
                                      a.noise_A_magnitude == b.noise_A_magnitude
