"""
`binsum` : this is a placeholder for a more general function, 
`binsum`, which still needs to be written!

The current version supports  `B = simple_binsum(i, A)`

* if `i` is a vector of integers and `A` is a vector, then `B` is a vector
of length `maximum(i)` with `B[n]` containing the sum of entries `A[t]` such that
`i[t] == n`.

* if `A` is a matrix, then the operation is applied column-wise, i.e., the
second dimension is summed over
"""
function simple_binsum{TI <: Integer, TF <: AbstractFloat}(i::Vector{TI},
                                                    A::Matrix{TF})
    if size(A, 2) != length(i)
        error("binsum: need size(A,2) = length(i)")
    end
    dim = size(A, 1)
    B = zeros(TF, dim, maximum(i))
    for m = 1:dim
        # @inbounds @simd
        for n = 1:length(i)
            B[m, i[n]] = B[m,i[n]] + A[m, n]
        end
    end
    return B
end

function simple_binsum{TI <: Integer, TF <: AbstractFloat}(i::Vector{TI},
                                                    A::Vector{TF})
    if length(A) != length(i)
        error("binsum: needs length(A) = length(i)")
    end
    B = zeros(TF, maximum(i))
    # this ought to be a SIMD loop. but that gives me a wrong answer! why?!
    for n = 1:length(i)
        B[i[n]] += + A[n]
    end
    return B
end

function binsum{TI <: Integer, TF <: AbstractFloat, N}(i::Vector{TI},
                                                       A::Array{TF, N},
                                                       dim::Integer)
    if length(i) != size(A, dim)
        error("binsum: length(i) must match size(A, dim)")
    end
    szA = size(A)
    nd_B = maximum(i)
    szB = [szA...]; szB[dim] = nd_B
    B = zeros(szB)
    # cartesian to linear indexing formula is
    #    i1 + (i2-1)*n1 + (i3-1) * n1*n2 + ... + (id-1) * n1*n2*...*n{d-1} + ...
    # so a linear index n can be corrected as follows:
    #    m = n1*n2*...*n{d-1}
    #    div(n, m) >>> (id-1) + (i{d+1}-1) * nd + ...
    #    mod( div(n, m), nd ) + 1 = id

    m = prod(szA[1:dim-1])
    nd = size(A, dim)

    for nA = 1:length(A)
        id = mod( div(nA, m), nd ) + 1
        nB = nA - m * (id-1) + m * (i[id]-1)
        B[nB] += A[nA]
    end
end


binsum{TI <: Integer, TF <: AbstractFloat}(i::Vector{TI},
                                           A::Vector{TF}) = binsum(i, A, 1)

colbinsum{TI <: Integer, TF <: AbstractFloat}(i::Vector{TI},
                                              A::Matrix{TF}) = binsum(i, A, 2)

rowbinsum{TI <: Integer, TF <: AbstractFloat}(i::Vector{TI},
                                              A::Matrix{TF}) = binsum(i, A, 1)
