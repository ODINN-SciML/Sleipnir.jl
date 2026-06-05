
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
        noise_A_magnitude::F = 5e-18
) where {F <: AbstractFloat}
    # Build PhysicalParameters based on values
    ft = typeof(g)
    physical_parameters = PhysicalParameters{ft}(ρ, g, ϵ, η₀,
        maxA, minA,
        maxC, minC,
        maxTlaw, minTlaw,
        noise_A_magnitude)

    return physical_parameters
end

function Base.:(==)(a::PhysicalParameters, b::PhysicalParameters)
    a.ρ == b.ρ && a.g == b.g &&
        a.ϵ == b.ϵ && a.η₀ == b.η₀ &&
        a.maxA == b.maxA && a.minA == b.minA &&
        a.minC == b.minC && a.maxC == b.maxC &&
        a.minTlaw == b.minTlaw && a.maxTlaw == b.maxTlaw &&
        a.noise_A_magnitude == b.noise_A_magnitude
end

# Show helpers
# Don't define them as closures over io, otherwise serialization with multiprocessing will fail
label(io, s, pad) = printstyled(io, rpad(s, pad); color = 183)
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
