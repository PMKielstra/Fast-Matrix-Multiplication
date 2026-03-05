module Common

function cannot_split(X, sz)
    any(size(X) .< sz)
end

function split_2d(X, sz)
    r_len, c_len = size(X)    
    r_step = r_len ÷ sz[1]
    c_step = c_len ÷ sz[2]
    r_ranges = [((i-1)*r_step + 1):(i == sz[1] ? r_len : i*r_step) for i in 1:sz[1]]
    c_ranges = [((j-1)*c_step + 1):(j == sz[2] ? c_len : j*c_step) for j in 1:sz[2]]
    [view(X, rr, cr) for rr in r_ranges, cr in c_ranges]
end

export make_recursive_multiplier

function make_recursive_multiplier(multiplier, A_shape, B_shape)
    function rmul(A, B)
        if cannot_split(A, A_shape) || cannot_split(B, B_shape)
            return A * B
        end
        next_A = split_2d(A, A_shape)
        next_B = split_2d(B, B_shape)
        multiplier(rmul, next_A, next_B)
    end
    rmul
end

export test_multiply

function test_multiply(multiplier, n1, n2)
    A = rand(1:100, n1, n2)
    B = rand(1:100, n2, n1)
    all(A * B == multiplier(A, B))
end

end