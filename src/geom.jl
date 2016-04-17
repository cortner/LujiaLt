
export Domain


#
# here we collect some functions that generate and manipulate
# simple lattice geometries, add defects, etc.
#

const Atri=[1.0 cos(pi/3); 0.0 sin(pi/3)]


"""Basic Lujia-Lt geometry type
* `X` : physical reference coordinates (both atomistic and FEM)
* `Z` : lattice index (integer) reference coordinates
* `mark`: 0 = atomistic, > 0 atomistic, but not core, < 0 continuum node;
     this gives users the flexibility to give more information to the markers
* `A` : lattice matrix

## Notes

*  X[:, 1:nA] = A * Z are the atomistic nodes; in particular note that FEM nodes
    don't have a corresponding Z entry
"""
type Domain
    X::Matrix{Float64}
    Z::Matrix{Int32}
    mark::Vector{Int8}
    A::Matrix{Float64}
    nA::Int
    tri::FEM.Triangulation
    info::Dict
end


# usual trick to allow more parameters to be stored
Base.getindex(geom::Domain, idx) = geom.info[idx]
Base.setindex!(geom::Domain, val, idx) = (geom.info[idx] = val)
   
"auxiliary function ala meshgrid"
2grid{T}(x::Vector{T}, y::Vector{T}) = x .+ ones(T, y)', ones(T, x) .+ y'

"turn 2grid into a list of points"
function 2grid_list(x, y)
    X, Y = 2grid(x, y)
    return [X[:]'; Y[:]']
end

"compute distance of points from some point"
dist(X, X0) = sqrt(sumabs2(X .- X0, 1));
dist(X) = sqrt(sumabs2(X, 1))


"""
generates a set of points X (2 x N) containing precisely (A Zᵈ) ∩ B(R).

## Parameters

All parameters are keyword. When `A != nothing` then `lattice` and `defect` are
ignored.

* `A` : lattice matrix
* `Ra` : {5.0} atomistic region radius; must be > 0
* `Rc` : {0.0} continuum region radius, if < Ra, then there is no continuum region
* `Rbuf` : {0.0} atomistic layer added to the atomistic core region B(Ra) >> B(Ra+Rbuf)
* `shape` : {:ball} this is the only shape allowed at the moment
* `lattice` : {:triangular}; this is the only admissible choice right now
* `defect` : {:none}; admissible choices are `:vacancy`, `:interstitial`
* `meshparams` : {[1.5; 3.0]}, 1.5 is the idea coarsening rate, 2.0 is the 
    maximum factor by which neighbouring layers of elements may increase

## Output

"""
function Domain(; A=nothing, Ra=5.0, Rc=0.0, Rbuf=0.0,
                lattice=:triangular, defect=:none, shape=:ball,
                meshparams = [1.5; 3.0] )

    if shape != :ball
        error("Domain: `shape` parameter allows only `:ball` at present")
    end
    
    # get the lattice matrix
    if A == nothing
        if lattice == :triangular
            AA = Atri
        else
            error("lattice_ball : only triangular lattice supported")
        end
    else
        AA = A
    end
    
    # generate a cubic portion of Zd containing -R/sig : R/sig in each
    # coordinate direction
    sig = minimum(svd(A)[2])
    ndim = ceil(Int, (Ra+Rbuf)/sig)
    geom.I = zeros(Int, 2*ndim+1, 2*ndim+1)
    Z = 2grid_list(-ndim:ndim, -ndim:ndim)
    # convert to physical reference configuration
    X = AA * Z
    # keep only the points inside Ra+Rbuf
    I = find(dist(X) .<= Ra+Rbuf)
    X = X[:, I]
    Z = Z[:, I]
    # collect some information
    nA = length(I)    
    Ia_core = find( dist(X) .<= Ra )
    mark = ones(int8, nA)    
    mark[Ia_core] = 0
    
    # ============================
    # now the continuum component
    # ============================
    if Rc > Ra+Rbuf
        R0 = Ra+Rbuf                           # initial radius
        H0 = minimum( setdiff(dist(X), 0.0) )  # initial mesh-size
        # create an array of radii and mesh sizes that we can adjust to the
        # requires Rc, before doing the actual assembly
        #####   <<<<<<<<  TODO
        #####
        r = R0
        h = H0
        layer = 0
        X = X[:]      # temporarily linearise the X array
        while r < Rc
            r += h
            layer -= 1
            # the ring has diameter 2 π r
            diam = 2 * π * r;
            # number of edges along this circle is ceil(diam / h)
            Nedges = ceil(diam/h)
            # create the new coordinates
            θ = linspace(0, diam - diam/Nedges, Nedges)
            Xnew = [r * cos(θ); r * sin(θ)]
            mark_new = layer * ones(Int8, length(θ))
            # append to the existing arrays
            append!(X, Xnew[:])
            append!(mark, mark_new)
        end
    end

    # create a triangulation
    X = reshape(X, 2, length(X)/2)
    tri = FEM.Triangulation(X)
    
    # create the reference domain
    geom = Domain(X, Z, mark, AA, nA, tri, Dict())

    # construct some defects (if requested by the user)
    if A == nothing
        # note that we are currently assuming that the lattice is
        # the triangulat lattice!
        if defect == :vacancy
            I0 = find( dist(X .< 0.1) )
            remove_atom!(geom, I0)
        elseif defect == :interstitial
            add_atom!(geom, [0.5;0.0])
        elseif defect != :none
            error("Domain : unknown `defect` parameter")
        end
    end

    return geom
end


function draw(geom::Domain; Y=nothing, U=nothing)
    
    
    
end


