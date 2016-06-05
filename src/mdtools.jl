

module MDTools

import Base.Cartesian.@nexprs

export NeighbourList


"""`type NeighbourList`

Implements a fast neighbourlist computation, loosely based on
[https://github.com/libAtoms/matscipy](MatSciPy).

### Usage

#### Basic Usage
```
X = rand(2, 10_000) * 100   # random positions in 100 x 100 cell
nlist = NeighbourList(X, 3.0)    # neig-list with cutoff 3.0
neigs = neighbours(nlist, 1)     # list of neighbours of atom X[:,1]
```

#### Iterators

*TODO*

#### Vectorised Usage
```
# assume X is a collection of points
nlist = NeighbourList(X, 3.0)
E = sum( nlist.r.^(-12) - 2.0 * nlist.r.^(-6) )
dE = zeros(size(X, 1), length(nlist))
dE[:, nlist.i] -= (-12) * nlist
# TODO: need a good binsum implementation
```

"""
type NeighbourList
    i::Vector{Int}
    j::Vector{Int}
    R::Matrix{Float64}    # this need not be Float64!!!!!
    r::Vector{Float64}    # but this probably does
    first::Vector{Int}
end

@inline binindex(n, ncells::NTuple{2, Int}, X, minX, L) =
    ( floor(Int, ncells[1] * (X[1, n] - minX[1])/L[1]) + 1,
      floor(Int, ncells[2] * (X[2, n] - minX[2])/L[2]) + 1 )

# @inline binindex{DIM}(n, ncells::NTuple{DIM, Int}, X, minX, L) =
#     ntuple(i->floor(Int, ncells[i]*(X[i,n]-minX[i])/L[i]), DIM)
#     # ( floor(Int, ncells[1] * (X[1, n] - minX[1])/L[1]) + 1,
#     #   floor(Int, ncells[2] * (X[2, n] - minX[2])/L[2]) + 1 )

@inline binindex_linear{DIM}(n, ncells::NTuple{DIM, Int}, X, minX, L) =
    sub2ind(ncells, binindex(n, ncells, X, minX, L)...)

# this should probably become a generated function!
@inline function neig_bins!(ncells::NTuple{2, Int}, bin::Int,
                            neig_bins::Vector{Int})
    nc1 = ncells[1]
    for a = -1:1
        neig_bins[2+a] = bin-nc1+a
        neig_bins[5+a] = bin+a
        neig_bins[8+a] = bin+nc1+a
    end
    return neig_bins
end

"returns number of (i,j) pairs stored in the neighbourlist"
Base.length(nlist::NeighbourList) = length(nlist.i)

"return number of atoms"
nX(nlist::NeighbourList) = length(nlist.first) - 1


function NeighbourList(X::Matrix{Float64}, cutoff::Float64)
    # analyse the configuration : compute extents of cluster,
    # then calculate how many bins are needed
    dim = size(X, 1)
    nX = size(X, 2)
    minX = zeros(dim)
    maxX = zeros(dim)
    L = zeros(dim)
    ncells = zeros(Int, dim)
    for i = 1:dim
        minX[i], maxX[i] = extrema(sub(X, i, :))
        # the 1.0001 ensures that all atoms are strictly inside the upper face
        L[i] = 1.0001 * (maxX[i] - minX[i])
        ncells[i] = max(1, floor(Int, L[i] / cutoff))
    end
    ncells_tup = tuple(ncells...)

    # now call an inner assembly function, which allows us to dispatch
    # on the type of ncells_tup !!!! this way we can make the dimension
    # a function parameter
    return assemble_nlist_inner(dim, nX, minX, L, ncells_tup, X, cutoff)
end


function assemble_nlist_inner{DIM}(dim::Int, nX::Int,
                                   minX::Vector{Float64},
                                   L::Vector{Float64},
                                   ncells::NTuple{DIM, Int},
                                   X::Matrix{Float64}, cutoff::Float64)

    ncells_tot = prod(ncells)

    # allocate space for a CCS-type linked-list
    seed = zeros(Int, ncells)   # the first atom in each bin
    last = zeros(Int, ncells)   # the currently last atom in each bin
    next = zeros(Int, nX)       # linking them up ...

    # loop over all atoms and put them into bins
    for n = 1:nX
        # compute bin index from current position X[:, n]
        bin = binindex_linear(n, ncells, X, minX, L)

        # check whether it is a seed
        if seed[bin] == 0
            seed[bin] = n
            last[bin] = n
        else
            next[last[bin]] = n
            last[bin] = n
        end
    end

    # now loop over atoms again to determine the neighbours
    alloc_fac = 12
    rsq = zeros(alloc_fac*nX)     # store square or distances
    R = zeros(alloc_fac*DIM*nX)   # store vectorial distance vectors of paris
    i = zeros(Int, alloc_fac*nX)
    j = zeros(Int, alloc_fac*nX)
    first = zeros(Int, nX+1)
    cutoff_sq = cutoff^2
    npairs = 0
    neig_bins = zeros(Int, 3^DIM)
    R_nm = zeros(DIM)
    nresize = 0
    gc_enable(false)
    for n = 1:nX
        # get cell-index again
        bin = binindex_linear(n, ncells, X, minX, L)

        for nbin in neig_bins!(ncells, bin, neig_bins)
            # skip this bin if it is outside the cell
            if nbin < 1 || nbin > ncells_tot; continue; end

            # loop through atoms in the nbin bin
            m = seed[nbin]
            @inbounds while m > 0
                if m != n
                    r_nm2 = 0.0
                    # need to use a generated function for this, probably not
                    # worthwhile?
                    # @nexprs DIM a->( R_nm[a] = X[a,m]-X[a,n];
                    #                  r_nm2 += R_nm[a]*R_nm[a] )
                    for a = 1:DIM   # should probably be unrolled
                        R_nm[a] = X[a,m]-X[a,n]
                        r_nm2 += R_nm[a]*R_nm[a]
                    end
                    if r_nm2 < cutoff_sq
                        npairs += 1
                        # check if arrays need to be resized
                        if npairs > length(i)
                            nresize += 1
                            resize!(i, 2*npairs)
                            resize!(j, 2*npairs)
                            resize!(rsq, 2*npairs)
                            resize!(R, 2*2*npairs)
                        end
                        # add new pair
                        i[npairs] = n
                        j[npairs] = m
                        rsq[npairs] = r_nm2
                        a0 = 2*npairs-2
                        # @nexprs DIM a->( R[a0+a] = R_nm[a] )
                        for a = 1:DIM
                            R[a0+a] = R_nm[a]
                        end
                        # remember the index for the first neighbour
                        # pair involving the centre-atom `n`
                        # this can be used to loop over neighbours
                        if first[n] == 0
                            first[n] = npairs
                        end
                    end
                end
                # next atom in the current neighbour-bin
                m = next[m]
            end  # end loop over atoms in neighbouring bin
        end # loop over bins
    end # loop over atoms

    # fix the `first` array in case one of the atoms has no neigbours at all?!
    # this makes sure that, if n has no neighbours, then first[n] == first[n+1]
    # and as a result the range first[n]:first[n+1] is empty.
    lastfirst = npairs
    for n = nX+1:-1:1
        if first[n] == 0
            first[n] = lastfirst
        else
            lastfirst = first[n]
        end
    end

    gc_enable(true)

    # postprocess the arrays
    # Rm = reshape(R, dim, length(i))
    Rm = reshape(R, DIM, length(i))[:, 1:npairs]
    i = i[1:npairs]
    j = j[1:npairs]
    rsq = rsq[1:npairs]

    # return the completed neighbour-list
    return NeighbourList(i, j, Rm, sqrt(rsq), first)
end


# #################### access functions for the neighbourlist  #################

export neighbours
"""basic access function to list of neighbours; better to iterature
via `sites(nlist)`; probably remove this!
"""
neighbours(nlist, n) = nlist.j[nlist.first[n]:nlist.first[n+1]-1]


export sites
# a simpler constructor, directly from an atoms object
sites(nlist::NeighbourList) = nlist

Base.start(nlist::NeighbourList) = 0::Int
Base.done(nlist::NeighbourList, state::Int) = (state == nX(nlist))

function Base.next(nlist::NeighbourList, state::Int)
    state += 1
    inds = nlist.first[state]:nlist.first[state+1]-1
    ret_tuple = (state,
                 nlist.j[inds],
                 nlist.r[inds],
                 nlist.R[:, inds]
                 )
    return ret_tuple, state
end


function uniquepairs(nlist::NeighbourList)
   ij = sort([nlist.i nlist.j], 2)
   ij = unique(ij, 1)
   return (ij[:,1], ij[:, 2])
end


end
