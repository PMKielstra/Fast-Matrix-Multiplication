"Convert a factorization of the matrix multiplication tensor ``<k, m, n>`` into an algorithm for multiplying matrices."
module FactorBased

"For a stack of matrices A, and a matrix B, compute a batched dot product.  In practice B will be a matrix of matrices itself, but we don't need to know that."
function two_by_four_dot(A, B)
    sz_A = size(A)
    sz_B = size(B)
    A_reshaped = reshape(A, sz_A[1], sz_A[2] * sz_A[3])
    B_reshaped = vec(B)
    function nonzero_sum(coeffs)
        nonzeros = (coeffs .!= 0)
        nonzero_coeffs = coeffs[nonzeros]
        sum(nonzero_coeffs .* B_reshaped[nonzeros])
    end
    reshape(nonzero_sum.(eachrow(A_reshaped)), sz_A[1])
end

export matmul_from_factorization

"""Given a list of arrays, returns a list of sub-arrays such that every slice in the first dimension of every array has some nonzero element."""
function filter_nonzero(arrays)
    function nonzero_indices(array)
        any.(eachslice(array .!= 0, dims=1))
    end
    indices = foldl((p, q) -> p .&& q, nonzero_indices.(arrays))
    (p -> p[indices, :, :]).(arrays)
end

"""
Builds a function to multiply a ``k \times m`` matrix and a ``m \times n`` one, given that ``<k, m, n> = \\sum_i \\alpha_i \\otimes \\beta_i \\otimes \\gamma_i``.
"""
function matmul_from_factorization(alphas, betas, gammas)
    alphas, betas, gammas = filter_nonzero((alphas, betas, gammas))
    function multiply(mul, A, B)
        dotted_As = two_by_four_dot(alphas, A)
        dotted_Bs = two_by_four_dot(betas, B)
        sum(kron(g, mul(dA, dB)) for (dA, dB, g) in zip(dotted_As, dotted_Bs, eachslice(gammas, dims=1)))
    end
end

end