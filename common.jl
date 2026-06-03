module Common

function cannot_split(X, sz)
    any(size(X) .< sz)
end

function split_2d(X, sz)
    r_len, c_len = size(X)
    r_step = r_len ÷ sz[1]
    c_step = c_len ÷ sz[2]
    r_ranges = [((i-1)*r_step+1):(i==sz[1] ? r_len : i*r_step) for i = 1:sz[1]]
    c_ranges = [((j-1)*c_step+1):(j==sz[2] ? c_len : j*c_step) for j = 1:sz[2]]
    [view(X, rr, cr) for rr in r_ranges, cr in c_ranges]
end

"""Pad matrix with zeros on the bottom and right so that it divides evenly into blocks of the given shape."""
function pad_to_divisibility(matrix, shape)
    mat_shape = size(matrix)
    pad_rows = mod(-mat_shape[1], shape[1]) # Smallest amount to pad such that the dimension is a multiple of 
    pad_cols = mod(-mat_shape[2], shape[2])

    if pad_rows == 0 && pad_cols == 0
        padded = matrix
    else
        padded = [matrix zeros(Int, mat_shape[1], pad_cols); zeros(Int, pad_rows, mat_shape[2]) zeros(Int, pad_rows, pad_cols)]
    end
    padded, pad_rows, pad_cols
end

export make_recursive_multiplier

function make_recursive_multiplier(multiplier, m, n, k)
    A_shape = (m, n)
    B_shape = (n, k)
    function rmul(A, B)
        if cannot_split(A, A_shape) || cannot_split(B, B_shape)
            return A * B
        end
        padded_A, remove_x, _ = pad_to_divisibility(A, A_shape)
        padded_B, _, remove_y = pad_to_divisibility(B, B_shape)
        next_A = split_2d(padded_A, A_shape)
        next_B = split_2d(padded_B, B_shape)
        result = multiplier(rmul, next_A, next_B)
        result[1:end-remove_x, 1:end-remove_y]
    end
end

export test_multiply

function test_multiply(multiplier, n)
    A = rand(1:100, n, n)
    B = rand(1:100, n, n)
    all(A * B == multiplier(A, B))
end

end
