
"""
A structure representing physical parameters used in simulations.

    PhysicalParameters{F <: AbstractFloat}

# Fields
- `ρ::F`: Density of ice.
- `g::F`: Gravitational acceleration.
- `ϵ::F`: A small parameter, often used for perturbations.
- `η₀::F`: Initial viscosity.
- `maxA::F`: Maximum A.
- `minA::F`: Minimum A.
- `maxC::F`: Maximum C.
- `minC::F`: Minimum C.
- `maxTlaw::F`: Maximum temperature according to some law.
- `minTlaw::F`: Minimum temperature according to some law.
- `noise_A_magnitude::F`: Magnitude of noise in A.
- `topo_spatial_window::F`: Spatial window for processing topography in meters, used for computing roughness and other inputs.
"""
struct PhysicalParameters{F <: AbstractFloat, I <: Integer} <: AbstractParameters
    ρ::F
    g::F
    ϵ::F
    η₀::F
    maxA::F
    minA::F
    maxC::F
    minC::F
    maxTlaw::F
    minTlaw::F
    noise_A_magnitude::F
    topo_spatial_window::F
    PDD_time_window::I
end

"""
Initialize the physical parameters of a model.

    PhysicalParameters(;
        ρ::Float64 = 900.0,
        g::Float64 = 9.81,
        ϵ::Float64 = 1e-3,
        η₀::Float64 = 1.0, 
        maxA::Float64 = 8e-17,
        minA::Float64 = 8.5e-20,
        maxC::Float64 = 1.0,
        minC::Float64 = 0.0,
        maxTlaw::Float64 = 1.0,
        minTlaw::Float64 = -25.0,
        noise_A_magnitude::Float64 = 5e-18,
        topo_spatial_window::Float64 = 200.0,
        PDD_time_window::Int64 = 7
        )

Keyword arguments
=================
    - `ρ`: Ice density
    - `g`: Gravitational constant
    - `ϵ`: Small number
    - `η₀`:  
    - `maxA`: Maximum value for `A` (Glen's coefficient)
    - `minA`: Minimum value for `A` (Glen's coefficient)
    - `maxC`: Maximum value of sliding coefficient `C`
    - `minC`: Minimum value of sliding coefficient `C`
    - `maxTlaw`: Maximum value of Temperature used in simulations on fake law
    - `minTlaw`: Minimum value of Temperature used in simulations on fake law
    - `noise_A_magnitude`: Magnitude of noise added to A
    - `topo_spatial_window`: Spatial window for processing topography in meters, used for computing roughness and other inputs.
    - `PDD_time_window`: Time window for positive degree days (PDD) in days, used to compute the CPDDs for inputs of laws.
"""
function PhysicalParameters(;
            ρ::F = 900.0,
            g::F = 9.81,
            ϵ::F = 1e-3,
            η₀::F = 1.0, 
            maxA::F = 8e-17,
            minA::F = 8.5e-20,
            maxC::F = 1.0,
            minC::F = 0.0,
            maxTlaw::F = 1.0,
            minTlaw::F = -25.0,
            noise_A_magnitude::F = 5e-18,
            topo_spatial_window::F = 200.0,
            PDD_time_window::I = 7
            ) where {F <: AbstractFloat, I <: Integer}
    # Build PhysicalParameters based on values
    ft = typeof(g)
    it = typeof(PDD_time_window)
    physical_parameters = PhysicalParameters{ft, it}(ρ, g, ϵ, η₀,
                                            maxA, minA,
                                            maxC, minC,
                                            maxTlaw, minTlaw,
                                            noise_A_magnitude,
                                            topo_spatial_window,
                                            PDD_time_window)

    return physical_parameters
end

Base.:(==)(a::PhysicalParameters, b::PhysicalParameters) = a.ρ == b.ρ && a.g == b.g && 
                                      a.ϵ == b.ϵ && a.η₀ == b.η₀ &&
                                      a.maxA == b.maxA && a.minA == b.minA && 
                                      a.minC == b.minC && a.maxC == b.maxC && 
                                      a.maxTlaw == b.maxTlaw &&
                                      a.noise_A_magnitude == b.noise_A_magnitude && a.topo_spatial_window == b.topo_spatial_window && 
                                      a.PDD_time_window == b.PDD_time_window 