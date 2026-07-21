
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
  - `temp_bias_min::F`: Minimum temperature bias (°C) for TI model calibration. Default: -10.0.
  - `temp_bias_max::F`: Maximum temperature bias (°C) for TI model calibration. Default: 10.0.
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
    temp_bias_min::F
    temp_bias_max::F
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
    - `temp_bias_min`: Minimum temperature bias (°C) for TI model calibration. Default: -10.0.
    - `temp_bias_max`: Maximum temperature bias (°C) for TI model calibration. Default: 10.0.
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
        prcp_fac_max::F = 10.0,
        temp_bias_min::F = -10.0,
        temp_bias_max::F = 10.0
) where {F <: AbstractFloat}
    ft = typeof(g)
    physical_parameters = PhysicalParameters{ft}(ρ, g, ϵ, η₀,
        maxA, minA,
        maxC, minC,
        maxTlaw, minTlaw,
        noise_A_magnitude,
        ρ_w, DDF_min, DDF_max, prcp_fac_min, prcp_fac_max,
        temp_bias_min, temp_bias_max)

    return physical_parameters
end

function Base.:(==)(a::PhysicalParameters, b::PhysicalParameters)
    a.ρ == b.ρ && a.g == b.g &&
        a.ϵ == b.ϵ && a.η₀ == b.η₀ &&
        a.maxA == b.maxA && a.minA == b.minA &&
        a.minC == b.minC && a.maxC == b.maxC &&
        a.minTlaw == b.minTlaw && a.maxTlaw == b.maxTlaw &&
        a.noise_A_magnitude == b.noise_A_magnitude &&
        a.ρ_w == b.ρ_w &&
        a.DDF_min == b.DDF_min && a.DDF_max == b.DDF_max &&
        a.prcp_fac_min == b.prcp_fac_min && a.prcp_fac_max == b.prcp_fac_max &&
        a.temp_bias_min == b.temp_bias_min && a.temp_bias_max == b.temp_bias_max
end

# Show helpers
# Don't define them as closures over io, otherwise serialization with multiprocessing will fail
label(io, s, pad) = printstyled(io, rpad(s, pad); color = 57)
sep(io) = printstyled(io, " · "; color = :light_black)
field(io, s) = printstyled(io, s; color = :light_black)
val(io, s) = print(io, s)
hint(io, s) = printstyled(io, s; color = :light_black)
check(b) = b ? "\e[32m✓\e[0m " : "\e[31m✗\e[0m "
nullable(io, x) = isnothing(x) ? hint(io, "(nothing)") : val(io, "$(nameof(typeof(x)))")

# Display setup
Base.show(io::IO, ::MIME"text/plain", params::PhysicalParameters) = Base.show(io, params)
function Base.show(io::IO, params::PhysicalParameters)
    pad = 12

    println(io, "PhysicalParameters")

    # Constants
    label(io, "  Constants", pad)
    field(io, "ρ");
    print(io, " = ");
    val(io, "$(params.ρ)");
    hint(io, " kg m⁻³")
    sep(io)
    field(io, "g");
    print(io, " = ");
    val(io, "$(params.g)");
    hint(io, " m s⁻²")
    sep(io)
    field(io, "η₀");
    print(io, " = ");
    val(io, "$(params.η₀)")
    println(io)

    # Glen A
    label(io, "  Glen A", pad)
    field(io, "min");
    print(io, " = ");
    val(io, "$(params.minA)")
    sep(io)
    field(io, "max");
    print(io, " = ");
    val(io, "$(params.maxA)");
    hint(io, " Pa⁻³ s⁻¹")
    println(io)

    # Sliding C
    label(io, "  Sliding C", pad)
    field(io, "min");
    print(io, " = ");
    val(io, "$(params.minC)")
    sep(io)
    field(io, "max");
    print(io, " = ");
    val(io, "$(params.maxC)")
    println(io)

    # Temp T
    label(io, "  Temp T", pad)
    field(io, "min");
    print(io, " = ");
    val(io, "$(params.minTlaw)");
    hint(io, " °C")
    sep(io)
    field(io, "max");
    print(io, " = ");
    val(io, "$(params.maxTlaw)");
    hint(io, " °C")
    println(io)

    # Numerics
    label(io, "  Numerics", pad)
    field(io, "ϵ");
    print(io, " = ");
    val(io, "$(params.ϵ)")
    println(io)
end
