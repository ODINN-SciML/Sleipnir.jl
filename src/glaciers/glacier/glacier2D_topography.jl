export compute_surface_topography, compute_surface_slope, compute_surface_aspect
export is_in_glacier

###############################################################
########  SURFACE TOPOGRAPHY & GRID MATH  #####################
###############################################################

function _window_cell_count(window_m::AbstractFloat, spacing::AbstractFloat)
    spacing_safe = max(spacing, Sleipnir.Float(eps(Float64)))
    n_cells = max(1, round(Int, window_m / spacing_safe))
    isodd(n_cells) ? n_cells : n_cells + 1
end

function _smoothed_surface(
        S::Matrix{<: AbstractFloat},
        Δx::AbstractFloat,
        Δy::AbstractFloat;
        window_m::AbstractFloat = 200.0)
    wx = _window_cell_count(window_m, Δx)
    wy = _window_cell_count(window_m, Δy)
    half_wx = (wx - 1) ÷ 2
    half_wy = (wy - 1) ÷ 2
    S_smooth = similar(S, Sleipnir.Float)

    for j in axes(S, 2), i in axes(S, 1)

        i0 = max(1, i - half_wx)
        i1 = min(size(S, 1), i + half_wx)
        j0 = max(1, j - half_wy)
        j1 = min(size(S, 2), j + half_wy)
        S_smooth[i, j] = Sleipnir.Float(mean(@view S[i0:i1, j0:j1]))
    end

    return S_smooth
end

function _centered_gradients(
        S::Matrix{<: AbstractFloat},
        Δx::AbstractFloat,
        Δy::AbstractFloat)
    nx, ny = size(S)
    dSdx = zeros(Sleipnir.Float, nx, ny)
    dSdy = zeros(Sleipnir.Float, nx, ny)

    if nx > 1
        dSdx_edges = Sleipnir.Float.(diff(S; dims = 1) ./ Δx)
        dSdx[1, :] .= dSdx_edges[1, :]
        dSdx[end, :] .= dSdx_edges[end, :]
        if nx > 2
            dSdx[2:(end - 1),
            :] .= 0.5 .*
                                    (dSdx_edges[1:(end - 1), :] .+ dSdx_edges[2:end, :])
        end
    end

    if ny > 1
        dSdy_edges = Sleipnir.Float.(diff(S; dims = 2) ./ Δy)
        dSdy[:, 1] .= dSdy_edges[:, 1]
        dSdy[:, end] .= dSdy_edges[:, end]
        if ny > 2
            dSdy[:,
            2:(end - 1)] .= 0.5 .*
                                    (dSdy_edges[:, 1:(end - 1)] .+ dSdy_edges[:, 2:end])
        end
    end

    return dSdx, dSdy
end

"""
    compute_surface_topography(
        S::Matrix{<: AbstractFloat},
        Δx::AbstractFloat,
        Δy::AbstractFloat;
        window_m::AbstractFloat = 200.0,
    )

Compute dynamic slope and aspect fields from the current glacier surface `S`.
The surface is spatially smoothed with a square moving window (`window_m`) before
computing finite-difference gradients to reduce pixel-scale noise.
"""
function compute_surface_topography(
        S::Matrix{<: AbstractFloat},
        Δx::AbstractFloat,
        Δy::AbstractFloat;
        window_m::AbstractFloat = 200.0)
    S_smooth = _smoothed_surface(S, Δx, Δy; window_m = window_m)
    dSdx, dSdy = _centered_gradients(S_smooth, Δx, Δy)

    slope = Sleipnir.Float.(rad2deg.(atan.(sqrt.(dSdx .^ 2 .+ dSdy .^ 2))))
    aspect = Sleipnir.Float.(mod.(rad2deg.(atan.(dSdy, -dSdx)), 360))
    return slope, aspect
end

function compute_surface_topography(
        glacier::Glacier2D;
        window_m::AbstractFloat = 200.0)
    return compute_surface_topography(
        glacier.S,
        glacier.Δx,
        glacier.Δy;
        window_m = window_m)
end

function compute_surface_slope(
        S::Matrix{<: AbstractFloat},
        Δx::AbstractFloat,
        Δy::AbstractFloat;
        window_m::AbstractFloat = 200.0)
    slope, _ = compute_surface_topography(S, Δx, Δy; window_m = window_m)
    return slope
end

function compute_surface_aspect(
        S::Matrix{<: AbstractFloat},
        Δx::AbstractFloat,
        Δy::AbstractFloat;
        window_m::AbstractFloat = 200.0)
    _, aspect = compute_surface_topography(S, Δx, Δy; window_m = window_m)
    return aspect
end

"""
    smooth!(A)

Smooths the interior of a 2D array `A` using a simple averaging method. The function modifies the array `A` in place.

# Arguments

  - `A::AbstractMatrix`: A 2D array to be smoothed.

# Details

The function updates the interior elements of `A` (excluding the boundary elements) by adding a weighted average of the second differences along both dimensions. The boundary elements are then set to the values of their nearest interior neighbors to maintain the boundary conditions.
"""
@views function smooth!(A)
    A[2:(end - 1),
    2:(end - 1)] .= A[2:(end - 1), 2:(end - 1)] .+
                                   1.0 ./ 4.1 .*
                                   (diff(diff(A[:, 2:(end - 1)], dims = 1), dims = 1) .+
                                    diff(diff(A[2:(end - 1), :], dims = 2), dims = 2))
    A[1, :]=A[2, :];
    A[end, :]=A[end - 1, :];
    A[:, 1]=A[:, 2];
    A[:, end]=A[:, end - 1]
end

# function smooth(A)
#     A_smooth = A[2:end-1,2:end-1] .+ 1.0./4.1.*(diff(diff(A[:,2:end-1], dims=1), dims=1) .+ diff(diff(A[2:end-1,:], dims=2), dims=2))
#     @tullio A_smooth_pad[i,j] := A_smooth[pad(i-1,1,1),pad(j-1,1,1)] # Fill borders 
#     return A_smooth_pad
# end

"""
    is_in_glacier(A::Matrix{F}, distance::I) where {I <: Integer, F <: AbstractFloat}

Return a matrix with booleans indicating if a given pixel is at distance at least
`distance` in the set of non zero values of the matrix. This usually allows
discarding the border pixels of a glacier.
A positive value of `distance` indicates a measurement from inside the glacier, while a negative `distance` indicates one from outside.

Arguments:

  - `A::Matrix{F}`: Matrix from which to compute the matrix of booleans.
  - `distance::I`: Distance to the border, computed as the number of pixels we need
    to move from within the glacier to find a pixel with value zero.
"""
function is_in_glacier(A::Matrix{F}, distance::I) where {I <: Integer, F <: AbstractFloat}
    B = convert.(F, (A .!= 0))
    # Reverse values in case we want distance from outside the border
    if distance < 0
        distance = -distance
        B .= 1.0 .- B
    end
    for i in 1:distance
        # We cannot use in-place affectation because this function is differentiated by Zygote in ODINN
        B = min.(
            B,
            circshift(B, (1, 0)),
            circshift(B, (-1, 0)),
            circshift(B, (0, 1)),
            circshift(B, (0, -1))
        )
    end
    B_bool = B .> 0.001
    if distance >= 0
        return B_bool
    else
        return .!B_bool
    end
end

"""
    block_average_pad_edge(mat::Matrix{F}, n::Int) where {F <: AbstractFloat}

Downsamples a matrix by averaging `n x n` blocks, using edge-replication padding
when the matrix dimensions are not divisible by `n`.
Edge padding replicates the last row/column values to expand the matrix so that both
dimensions are divisible by `n`.
Returns a matrix of averaged values with size `(ceil(Int, X/n), ceil(Int, Y/n))`.

Arguments

  - `mat::Matrix{F}`: Input 2D matrix.
  - `n::Int`: Block size for downsampling.
"""
function block_average_pad_edge(mat::Matrix{F}, n::Int) where {F <: AbstractFloat}
    X, Y = size(mat)
    new_X = ceil(Int, X / n) * n
    new_Y = ceil(Int, Y / n) * n

    # Create padded matrix filled with edge values
    padded = similar(mat, new_X, new_Y)

    # Fill original data
    padded[1:X, 1:Y] .= mat

    # Pad bottom rows with last row
    if new_X > X
        for i in (X + 1):new_X
            padded[i, 1:Y] .= mat[end, :]
        end
    end

    # Pad right columns with last column
    if new_Y > Y
        for j in (Y + 1):new_Y
            padded[1:X, j] .= mat[:, end]
        end
    end

    # Fill bottom-right corner (if both X and Y were padded)
    if new_X > X && new_Y > Y
        for i in (X + 1):new_X
            for j in (Y + 1):new_Y
                padded[i, j] = mat[end, end]
            end
        end
    end

    return block_average(padded, n)
end

"""
    block_average(mat::Matrix{F}, n::Int) where {F <: AbstractFloat}

Downsamples a matrix by averaging non-overlapping `n x n` blocks.
Returns a matrix of the block-averaged values with size `(div(X, n), div(Y, n))`
where `(X, Y) = size(mat)`.

Arguments

  - `mat::Matrix{F}`: Input 2D matrix.
  - `n::Int`: Block size for downsampling. Both matrix dimensions must be divisible by `n`.
"""
function block_average(mat::Matrix{F}, n::Int) where {F <: AbstractFloat}
    X, Y = size(mat)
    @assert X % n == 0 && Y % n == 0 "Matrix dimensions are $(size(mat)) but they are not divisible by n=$n"

    A, B = div(X, n), div(Y, n)
    reshaped = reshape(mat, n, A, n, B)
    permuted = permutedims(reshaped, (2, 4, 1, 3))  # (A, B, n, n)
    mean_blocks = mean(permuted, dims = (3, 4))
    return dropdims(mean_blocks, dims = (3, 4))
end
