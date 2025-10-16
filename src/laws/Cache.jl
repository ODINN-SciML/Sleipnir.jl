import Base.similar, Base.size, Base.fill
export Cache
export MatrixCacheNoVJP, ScalarCacheNoVJP, MatrixCache, ScalarCache

"""
    Cache

Abstract type for defining a cache struct of a law.

Mandatory field:
- `value`: Store the result of a forward evaluation.

Optional fields:
- `vjp_inp`: Store the result of the evaluation of the vector-Jacobian product (VJP) with respect to the inputs.
- `vjp_θ`: Store the result of the evaluation of the vector-Jacobian product (VJP) with respect to the parameters θ.

Notes:
- If a concrete subtype does not implement `vjp_inp` and `vjp_θ`, the law should not be used for gradient computation, and therefore for inversions.
"""
abstract type Cache end

"""
    MatrixCacheNoVJP <: Cache

A mutable cache structure for storing a two-dimensional array of `Float64` values.
This is typically used for spatially varying laws.
This struct is intended for use cases where the law is not differentiated, and hence the vector-Jacobian products (VJP) are not required.
Fields:
- `value::Array{Float64, 2}`: The cached matrix.
"""
mutable struct MatrixCacheNoVJP <: Cache
    value::Array{Float64, 2}
end

"""
    ScalarCacheNoVJP <: Cache

A mutable cache structure for storing a scalar value as a zero-dimensional array of `Float64`.
This is typically used for constant per glacier laws.
This struct is intended for use cases where the law is not differentiated, and hence the vector-Jacobian products (VJP) are not required.
Fields:
- `value::Array{Float64, 0}`: The cached scalar value.
"""
mutable struct ScalarCacheNoVJP <: Cache
    value::Array{Float64, 0}
end

similar(c::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}) = typeof(c)(similar(c.value))
size(c::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}) = size(c.value)
fill(c::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}, s) = typeof(c)(fill(c.value, s))
Base.:(==)(a::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}, b::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}) = a.value == b.value

"""
    MatrixCache <: Cache

A cache structure for storing a two-dimensional array of `Float64` values along with
their associated vector-Jacobian products (VJP).
This is typically used for spatially varying laws.
Fields:
- `value::Array{Float64, 2}`: The cached matrix.
- `vjp_inp::Array{Float64, 2}`: VJP with respect to inputs.
- `vjp_θ::Vector{Float64}`: VJP with respect to parameters.
"""
struct MatrixCache <: Cache
    value::Array{Float64, 2}
    vjp_inp::Array{Float64, 2}
    vjp_θ::Vector{Float64}
end

"""
    ScalarCache <: Cache

A cache structure for storing a scalar value as a zero-dimensional array of `Float64` along with
their associated vector-Jacobian products (VJP).
This is typically used for constant per glacier laws.
Fields:
- `value::Array{Float64, 0}`: The cached scalar value.
- `vjp_inp::Array{Float64, 0}`: VJP with respect to inputs.
- `vjp_θ::Vector{Float64}`: VJP with respect to parameters.
"""
struct ScalarCache <: Cache
    value::Array{Float64, 0}
    vjp_inp::Array{Float64, 0}
    vjp_θ::Vector{Float64}
end

similar(c::Union{ScalarCache, MatrixCache}) = typeof(c)(similar(c.value), similar(c.vjp_inp), similar(c.vjp_θ))
size(c::Union{ScalarCache, MatrixCache}) = size(c.value)
fill(c::Union{ScalarCache, MatrixCache}, s) = typeof(c)(fill(c.value, s), fill(c.vjp_inp, s), similar(c.vjp_θ))
Base.:(==)(a::Union{ScalarCache, MatrixCache}, b::Union{ScalarCache, MatrixCache}) = a.value == b.value && a.vjp_inp == b.vjp_inp && a.vjp_θ == b.vjp_θ
