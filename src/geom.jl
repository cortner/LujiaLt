
using Compose

#
# here we collect some functions that generate and manipulate
# simple lattice geometries, add defects, etc.
#

const Atri=[1.0 cos(pi/3); 0.0 sin(pi/3)]


# """Basic Lujia-Lt geometry type
# * `X` : physical reference coordinates (both atomistic and FEM)
# * `Z` : lattice index (integer) reference coordinates
# * `mark`: 0 = atomistic, > 0 atomistic, but not core, < 0 continuum node;
#      this gives users the flexibility to give more information to the markers
# * `A` : lattice matrix

# ## Notes

# *  X[:, 1:nA] = A * Z are the atomistic nodes; in particular note that FEM nodes
#     don't have a corresponding Z entry
# """
# type Domain
#     X::Matrix{Float64}
#     Z::Matrix{Int32}
#     mark::Vector{Int8}
#     A::Matrix{Float64}
#     nA::Int
#     tri::FEM.Triangulation
#     info::Dict
# end


# usual trick to allow more parameters to be stored
Base.getindex(geom::Domain, idx) = geom.info[idx]
Base.setindex!(geom::Domain, val, idx) = (geom.info[idx] = val)
   
"auxiliary function ala meshgrid"
twogrid{T}(x::Vector{T}, y::Vector{T}) =
    (x .+ ones(T, length(y))', ones(T, length(x)) .+ y')
# twogrid(x, y) = twogrid(collect(x), collect(y))

"turn 2grid into a list of points"
function twogrid_list(x, y)
    X, Y = twogrid(collect(x), collect(y))
    return [X[:]'; Y[:]']
end

"compute distance of points from some point"
dist(X, X0) = sqrt(sumabs2(X .- X0, 1));
dist(X) = sqrt(sumabs2(X, 1))
export dist


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
    sig = minimum(svd(AA)[2])
    ndim = ceil(Int, (Ra+Rbuf)/sig)
    Z = twogrid_list(-ndim:ndim, -ndim:ndim)
    # convert to physical reference configuration
    X = AA * Z
    # keep only the points inside Ra+Rbuf
    I = find(dist(X) .<= Ra+Rbuf)
    X = X[:, I]
    Z = Z[:, I]
    # collect some information
    nA = length(I)    
    Ia_core = find( dist(X) .<= Ra )
    mark = ones(Int8, nA)   
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
        h = H0
        r = R0-0.4*h
        layer = 0
        X = X[:]      # temporarily linearise the X array
        while r < Rc
            r += h
            layer -= 1
            # the current ring has diameter 2 π r
            diam = 2 * π * r;
            # number of edges along this circle is ceil(diam / h)
            Nedges = ceil(diam/h)
            # create the new coordinates
            θ = collect(linspace(0, 2*π - 2*π/Nedges, Nedges))
            Xnew = [r * cos(θ)'; r * sin(θ)']
            mark_new = layer * ones(Int8, length(θ))
            # append to the existing arrays
            append!(X, Xnew[:])
            append!(mark, mark_new)
            # update h
            hfun = (r/R0)^meshparams[1]
            h = max( (1.5*h), min(meshparams[1]*h, hfun) )
        end
    end
    
    # create a triangulation
    # TODO: CURRENTLY WE DO THIS TWICE!!!!! (again after introducing the defect)
    X = reshape(X, 2, length(X) ÷ 2)
    tri = FEM.Triangulation(X)
    
    # create the reference domain
    geom = Domain(X, Z, mark, AA, nA, tri, Dict())

    # construct some defects (if requested by the user)
    if A == nothing
        # note that we are currently assuming that the lattice is
        # the triangular lattice!
        if defect == :vacancy
            I0 = find( dist(X) .< 0.1 )
            remove_atom!(geom, I0[1])
        elseif defect == :interstitial
            add_atom!(geom, [0.5;0.0])
        elseif defect != :none
            error("Domain : unknown `defect` parameter")
        end
    end
    
    return geom
end


"number of nodes (not free nodes!)"
nX(geom::Domain) = size(geom.X, 2)



"remove an atom from the geometry, usually to create a vacancy"
function remove_atom!(geom::Domain, idx)
    if geom.mark[idx] != 0
        error("remove_atom! : only allowed to remove core atoms")
    end
    # TODO: if there are continuum nodes, then they are also included
    # but not in Z!!!!
    Iall = setdiff(1:nX(geom), idx)
    geom.X = geom.X[:, Iall]
    geom.mark = geom.mark[Iall]
    Iat = setdiff(find(geom.mark .>= 0), idx)
    geom.Z = geom.Z[:, Iat]
    geom.nA = geom.nA - 1
    # redo the triangulation
    geom.tri = FEM.Triangulation(geom.X)
end
export remove_atom!



"""
visualise the domain and if passed deformation or displacement
"""
function compose(geom::Domain;
                 Y=nothing, U=nothing,
                 axis=:auto, lwidth=:auto)
    if axis == :auto
        axis = Plotting.autoaxis(geom.X)
    end
    if lwidth == :auto
        lwidth = 0.3
    end
    # plot core atomistic nodes
    ctx_atm = Plotting.compose_atoms(geom.X[:, find(geom.mark .== 0)];
                                     axis=axis, fillcolor="tomato")
    # plot buffer nodes
    ctx_buf = Plotting.compose_atoms(geom.X[:, find(geom.mark .> 0)];
                                     axis=axis,
                                     fillcolor="aliceblue",
                                     linecolor="darkblue" )
    # ctx_cb = Plotting.compose_atoms(geom.X[:, find(geom.mark .< 0)];
    #                                 axis=axis, fillcolor="aliceblue")
    # plot the finite element mesh (should we remove the inner triangles)?
    ctx_cb = Plotting.compose_elements(geom.tri.X, geom.tri.T;
                                       axis=axis, lwidth=lwidth)
    return Compose.compose( context(), ctx_atm, ctx_buf, ctx_cb )
end



