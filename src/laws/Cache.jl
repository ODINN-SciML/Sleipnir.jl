import Base.similar, Base.size, Base.fill
export MatrixCacheNoVJP, ScalarCacheNoVJP, MatrixCache, ScalarCache

"""
    MatrixCacheNoVJP

A mutable cache structure for storing a two-dimensional array of `Float64` values.
This is typically used for spatially varying laws.
This struct is intended for use cases where the law is not differentiated, and hence the vector-Jacobian products (VJP) are not required.
Fields:
- `value::Array{Float64, 2}`: The cached matrix.
"""
mutable struct MatrixCacheNoVJP
    value::Array{Float64, 2}
end

"""
    ScalarCacheNoVJP

A mutable cache structure for storing a scalar value as a zero-dimensional array of `Float64`.
This is typically used for constant per glacier laws.
This struct is intended for use cases where the law is not differentiated, and hence the vector-Jacobian products (VJP) are not required.
Fields:
- `value::Array{Float64, 0}`: The cached scalar value.
"""
mutable struct ScalarCacheNoVJP
    value::Array{Float64, 0}
end

similar(c::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}) = typeof(c)(similar(c.value))
size(c::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}) = size(c.value)
fill(c::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}, s) = typeof(c)(fill(c.value, s))
Base.:(==)(a::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}, b::Union{ScalarCacheNoVJP, MatrixCacheNoVJP}) = a.value == b.value

"""
    MatrixCache

A cache structure for storing a two-dimensional array of `Float64` values along with
their associated vector-Jacobian products (VJP).
This is typically used for spatially varying laws.
Fields:
- `value::Array{Float64, 2}`: The cached matrix.
- `vjp_inp::Array{Float64, 2}`: VJP with respect to inputs.
- `vjp_θ::Vector{Float64}`: VJP with respect to parameters.
"""
struct MatrixCache
    value::Array{Float64, 2}
    vjp_inp::Array{Float64, 2}
    vjp_θ::Vector{Float64}
end

"""
    ScalarCache

A cache structure for storing a scalar value as a zero-dimensional array of `Float64` along with
their associated vector-Jacobian products (VJP).
This is typically used for constant per glacier laws.
Fields:
- `value::Array{Float64, 0}`: The cached scalar value.
- `vjp_inp::Array{Float64, 0}`: VJP with respect to inputs.
- `vjp_θ::Vector{Float64}`: VJP with respect to parameters.
"""
struct ScalarCache
    value::Array{Float64, 0}
    vjp_inp::Array{Float64, 0}
    vjp_θ::Vector{Float64}
end

similar(c::Union{ScalarCache, MatrixCache}) = typeof(c)(similar(c.value), similar(c.vjp_inp), similar(c.vjp_θ))
size(c::Union{ScalarCache, MatrixCache}) = size(c.value)
fill(c::Union{ScalarCache, MatrixCache}, s) = typeof(c)(fill(c.value, s), fill(c.vjp_inp, s), similar(c.vjp_θ))
Base.:(==)(a::Union{ScalarCache, MatrixCache}, b::Union{ScalarCache, MatrixCache}) = a.value == b.value && a.vjp_inp == b.vjp_inp && a.vjp_θ == b.vjp_θ
