export safe_approx, is_in_glacier

# Function to override "≈" to handle nothing values
function safe_approx(a, b)
    if isnothing(a) && isnothing(b)
        return true
    elseif isnothing(a) || isnothing(b)
        return false
    else
        return a ≈ b
    end
end

"""
    is_in_glacier(A::Matrix{F}, distance::I) where {I <: Integer, F <: AbstractFloat}

Return a matrix with booleans indicating if a given pixel is at distance at least
`distance` in the set of non zero values of the matrix. This usually allows
discarding the border pixels of a glacier.

Arguments:
- `A::Matrix{F}`: Matrix from which to compute the matrix of booleans.
- `distance::I`: Distance to the border, computed as the number of pixels we need
    to move to find a pixel with value zero.
"""
function is_in_glacier(A::Matrix{F}, distance::I) where {I <: Integer, F <: AbstractFloat}
    B = convert.(F, (A .!= 0))
    for i in 1:distance
        B .= min.(B, circshift(B, (1,0)), circshift(B, (-1,0)), circshift(B, (0,1)), circshift(B, (0,-1)))
    end
    return B .> 0.001
end
