
using Compose
using PyCall

export nX, positions, dist, remove_atom!


#
# here we collect some functions that generate and manipulate
# simple lattice geometries, add defects, etc.
# the `Domain` type is already defined in LujiaLt.jl
#

const Atri=[1.0 cos(pi/3); 0.0 sin(pi/3)]

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


"return a lattice ball"
function lattice_ball(;R=5.0, A=Atri, x0=[0.0;0.0])
    # generate a cubic portion of Zd containing -R/sig : R/sig in each
    # coordinate direction
    sig = minimum(svd(A)[2])
    ndim = ceil(Int, (norm(x0)+R)/sig)
    Z = twogrid_list(-ndim:ndim, -ndim:ndim)
    # convert to physical reference configuration
    X = A * Z
    # keep only the points inside R
    I = find(dist(X, x0) .<= R)
    X = X[:, I]
    Z = Z[:, I]
    return X, Z
end


"construct FEM nodes on an annulus"
function radial_nodes(Rin, Rout;
                      hin = 1.0, hgrowth=1.5, growthbound=3.0,
                      x0=[0.0;0.0])
   # arrays to be returned
   X = Float64[]    # positions of nodes
   mark = Int[]     # markers (layers)

   # create an array of radii and mesh sizes that we can adjust to the
   # requires Rc, before doing the actual assembly
   #####   <<<<<<<<  TODO

   # initialise
   r = Rin-hin
   h = hin
   layer = 0
   # loop over layers
   while r < Rout
      r += h
      layer -= 1
      # the current ring has diameter 2 π r
      diam = 2 * π * r;
      # number of edges along this circle is ceil(diam / h)
      Nedges = ceil(diam/h)
      # create the new coordinates
      θ = collect(linspace(0, 2*π - 2*π/Nedges, Nedges))
      Xnew = [x0[1] + r * cos(θ)'; x0[2] + r * sin(θ)']
      mark_new = layer * ones(Int8, length(θ))
      # append to the existing arrays
      append!(X, Xnew[:])
      append!(mark, mark_new)
      # update h
      hfun = (r/Rin)^hgrowth
      h = max( (1.5*h), min(growthbound*h, hfun) )
   end

   return reshape(X, 2, length(X) ÷ 2), mark
end




function Domain(; A=nothing, Ra=5.0, Rc=0.0,
                lattice=:triangular, defect=:none, shape=:ball,
                meshparams = [1.5; 3.0], x0 = :auto )

   if shape != :ball
      error("Domain: `shape` parameter allows only `:ball` at present")
   end

   # get the lattice matrix
   if A == nothing
      lattice == :triangular ? AA = Atri :
         error("lattice_ball : only triangular lattice supported")
   else
      AA = A
   end

    # default defect core position
    if x0 == :auto
      x0 = [0.0; 0.0]
      if defect == :screw
         x0 = [0.5; sqrt(3)/4]
      end
      if defect == :interstitial
         x0 = [0.5; 0.0]
      end
   end


    X, Z = lattice_ball(A=AA, R=Ra, x0=x0)
    # collect some information
    nA = size(X, 2)
    Ia_core = find( dist(X, x0) .<= Ra )
    mark = ones(Int8, nA)
    mark[Ia_core] = 0

   #  ============================
   #  now the continuum component
   #  ============================
   if Rc > Ra
      X_c, mark_c = radial_nodes(Ra + 0.6, Rc,
                                 hgrowth=meshparams[1],
                                 growthbound=meshparams[2])
      # add to existing array
      X = [X X_c]
      mark = [mark; mark_c]
   end

    # create a Triangulation object without actual triangulation
    # the proper triangulation will be created after we add the defects
    X = reshape(X, 2, length(X) ÷ 2)
    tri = FEM.Triangulation(X, Matrix{Int}(), PyObject(nothing))

    # create the reference domain
    geom = Domain(Z, mark, AA, nA, tri, Dict())

    # construct some defects (if requested by the user)
    if A == nothing
        # note that we are currently assuming that the lattice is
        # the triangular lattice!
        # remove_atom! and add_atom! both redo the triangulation
        # so we don't need to create a new one!
        if defect == :vacancy
            I0 = find( dist(X) .< 0.1 )
            remove_atom!(geom, I0[1])
        elseif defect == :interstitial
            add_atom!(geom, [0.5;0.0])
        elseif defect == :none
            geom.tri = FEM.Triangulation(X)
        else
            error("Domain : unknown `defect` parameter")
        end
    end

    return geom
end


"number of nodes (not free nodes!)"
nX(geom::Domain) = size(geom.tri.X, 2)

"reference to position array, same as `X`"
positions(geom::Domain) = geom.tri.X

"""
domain dimension - currently this must always be 2; this function is
included to ensure readability of some parts of the code. (and for the future
allow an extension maybe?)
"""
ddim(geom::Domain) = 2


"remove an atom from the geometry, usually to create a vacancy"
function remove_atom!(geom::Domain, idx)
   if geom.mark[idx] != 0
      error("remove_atom! : only allowed to remove core atoms")
   end
   Iall = setdiff(1:nX(geom), idx)
   X = positions(geom)[:, Iall]
   geom.mark = geom.mark[Iall]
   # note, geom.Z is shorter than positions(geom) since continuum nodes are
   # not included
   Iat = setdiff(find(geom.mark .>= 0), idx)
   geom.Z = geom.Z[:, Iat]
   geom.nA = geom.nA - 1
   # redo the triangulation
   geom.tri = FEM.Triangulation(X)
   return geom
end


"remove an atom from the geometry, usually to create a vacancy"
function add_atom!(geom::Domain, x0::Vector)
   Iall = setdiff(1:nX(geom), idx)
   X = [positions(geom) x0]
   geom.mark = [geom.mark; 0]
   geom.nA = geom.nA + 1
   # redo the triangulation
   geom.tri = FEM.Triangulation(X)
   return geom
end
