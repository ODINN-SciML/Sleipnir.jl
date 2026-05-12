
"""
A structure representing physical parameters used in simulations.

    PhysicalParameters{F <: AbstractFloat}

# Fields

  - `ρ::F`: Density of ice.
  - `g::F`: Gravitational acceleration.
  - `ϵ::F`: Regularization used in the square root of norms for AD numerical stability.
  - `η₀::F`: Initial viscosity.
  - `maxA::F`: Maximum A.
  - `minA::F`: Minimum A.
  - `maxC::F`: Maximum C.
  - `minC::F`: Minimum C.
  - `maxTlaw::F`: Maximum temperature according to some law.
  - `minTlaw::F`: Minimum temperature according to some law.
  - `noise_A_magnitude::F`: Magnitude of noise in A.
  - `ρ_w::F`: Density of water (kg m⁻³), used for ice-to-water-equivalent conversions.
  - `DDF_min::F`: Minimum degree-day factor for TI model calibration (m w.e. °C⁻¹ d⁻¹).
  - `DDF_max::F`: Maximum degree-day factor for TI model calibration (m w.e. °C⁻¹ d⁻¹).
  - `prcp_fac_min::F`: Minimum precipitation correction factor for TI model calibration.
  - `prcp_fac_max::F`: Maximum precipitation correction factor for TI model calibration.
"""
struct PhysicalParameters{F <: AbstractFloat} <: AbstractParameters
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
    ρ_w::F
    DDF_min::F
    DDF_max::F
    prcp_fac_min::F
    prcp_fac_max::F
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
        maxC::Float64 = 8e-17, # TODO: to be revised
        minC::Float64 = 8.5e-20,
        maxTlaw::Float64 = 1.0,
        minTlaw::Float64 = -25.0,
        noise_A_magnitude::Float64 = 5e-18
        )

# Keyword arguments

    - `ρ`: Ice density
    - `g`: Gravitational acceleration.
    - `ϵ`: Regularization used in the square root of norms for AD numerical stability.
    - `η₀`: Factor to cap surface elevation differences with the upstream ice thickness to impose boundary condition in the iceflow equation
    - `maxA`: Maximum value for `A` (Glen's coefficient)
    - `minA`: Minimum value for `A` (Glen's coefficient)
    - `maxC`: Maximum value of sliding coefficient `C`
    - `minC`: Minimum value of sliding coefficient `C`
    - `maxTlaw`: Maximum value of Temperature used in simulations on fake law
    - `minTlaw`: Minimum value of Temperature used in simulations on fake law
    - `noise_A_magnitude`: Magnitude of noise added to A
    - `ρ_w`: Water density (kg m⁻³). Default: 1000.0.
    - `DDF_min`: Minimum degree-day factor for TI model calibration (m w.e. °C⁻¹ d⁻¹). Default: 0.5×10⁻³.
    - `DDF_max`: Maximum degree-day factor for TI model calibration (m w.e. °C⁻¹ d⁻¹). Default: 20.0×10⁻³.
    - `prcp_fac_min`: Minimum precipitation correction factor for TI model calibration. Default: 0.1.
    - `prcp_fac_max`: Maximum precipitation correction factor for TI model calibration. Default: 10.0.
"""
function PhysicalParameters(;
        ρ::F = 900.0,
        g::F = 9.81,
        ϵ::F = 1e-10,
        η₀::F = 1.0,
        maxA::F = 8e-17,
        minA::F = 8.5e-20,
        maxC::F = 8e-17, # TODO: to be revised
        minC::F = 8.5e-20,
        maxTlaw::F = 1.0,
        minTlaw::F = -25.0,
        noise_A_magnitude::F = 5e-18,
        ρ_w::F = 1000.0,
        DDF_min::F = 0.5 / 1000.0,
        DDF_max::F = 20.0 / 1000.0,
        prcp_fac_min::F = 0.1,
        prcp_fac_max::F = 10.0
) where {F <: AbstractFloat}
    # Build PhysicalParameters based on values
    ft = typeof(g)
    physical_parameters = PhysicalParameters{ft}(ρ, g, ϵ, η₀,
        maxA, minA,
        maxC, minC,
        maxTlaw, minTlaw,
        noise_A_magnitude,
        ρ_w, DDF_min, DDF_max, prcp_fac_min, prcp_fac_max)

    return physical_parameters
end

function Base.:(==)(a::PhysicalParameters, b::PhysicalParameters)
    a.ρ == b.ρ && a.g == b.g &&
        a.ϵ == b.ϵ && a.η₀ == b.η₀ &&
        a.maxA == b.maxA && a.minA == b.minA &&
        a.minC == b.minC && a.maxC == b.maxC &&
        a.maxTlaw == b.maxTlaw &&
        a.noise_A_magnitude == b.noise_A_magnitude &&
        a.ρ_w == b.ρ_w &&
        a.DDF_min == b.DDF_min && a.DDF_max == b.DDF_max &&
        a.prcp_fac_min == b.prcp_fac_min && a.prcp_fac_max == b.prcp_fac_max
end
