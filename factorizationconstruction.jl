module FactorizationConstruction

"""List partitions of n of length l"""
function partitions_of_length(n, l)
    if l == 1
        return [[n]]
    end
    vcat([vcat([i], p) for i in 0:n for p in partitions_of_length(n - i, l - 1)])
end

"""Compute the Kronecker product of two matrices with entries in K[eps], written as a list of matrices in K."""
function power_series_kron(p1, p2)
    power1 = length(p1) - 1 # First element of p1 is eps^0
    power2 = length(p2) - 1
    power = power1 + power2 # Later elements are discarded as O(eps^(power + 1))
    [
        sum(
            kron(p1[i + 1], p2[current_power - i + 1]) for i in max(0, current_power - power2):min(power1, current_power), dims=1
        ) for current_power in 0:power
    ]
end

export square_up

"""Given a factorization of <m, n, k> in power series from K[eps], produce a factorization of <mnk, mnk, mnk>."""
function square_up(alphas, betas, gammas)
    pslk(l1, l2, l3) = [power_series_kron(power_series_kron(x[1], y[1]), z[1]) for x in eachslice(l1, dims=1) for y in eachslice(l2, dims=1) for z in eachslice(l3, dims=1)]
    pslk(alphas, betas, gammas), pslk(betas, gammas, alphas), pslk(gammas, alphas, betas)
end

export to_power

"""Given a factorization of <N, N, N> in power series from K[eps], produce a factorization of <N^s, N^s, N^s>."""
function to_power(s, alphas, betas, gammas)
    function _to_power_single(q)
        result = q
        for _ in 1:s
            result = power_series_kron(result, q)
        end
        result
    end
    to_power_single(alphas), to_power_single(betas), to_power_single(gammas)
end

export border_factorization_to_full_factorization

"""Given a vector `factors` of power series from K[eps], return a vector of factorizations with elements in K."""
function border_factorization_to_full_factorization(factors)
    l = length(factors)
    original_rk = length(factors[1])
    power = length(factors[1][1])
    partitions = partitions_of_length(power - 1, l)
    [
        [factors[i][r][p[i] + 1] for r in 1:original_rk for p in partitions] for i in 1:l
    ]
end

end