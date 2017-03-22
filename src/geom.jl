
# this is part of the main LujiaLt module

export nX, positions, dist, remove_atom!


# here we collect some functions that generate and manipulate
# simple lattice geometries, add defects, etc.
# the `Domain` type is already defined in LujiaLt.jl

const Atri=[1.0 cos(pi/3); 0.0 sin(pi/3)]

# usual trick to allow more parameters to be stored
Base.getindex(geom::Domain, idx) = geom.info[idx]
Base.setindex!(geom::Domain, val, idx) = (geom.info[idx] = val)

"auxiliary function ala meshgrid"
meshgrid{T}(x::Vector{T}, y::Vector{T}) =
   (x .+ ones(T, length(y))', ones(T, length(x)) .+ y')

"turn 2grid into a list of points"
function meshgrid_list(x, y)
   X, Y = meshgrid(collect(x), collect(y))
   return [X[:]'; Y[:]']
end

# "`dgrid(x, dim)` : another meshgrid like thing to replace meshgrid"
# function dgrid{T}(x::Vector{T}, dim)
#    o = ones(T, length(x))
#    dim == 2 && return x .+ o', o .* x'
#    dim == 3 && return kron(x, o, o), kron(o, x, o), kron(o, o, x)
#    throw(ArgumentError("dgrid: dim must be 2, or 3."))
# end
#
# "dgrid in a list of vectors form"
# dgrid_list(x, dim) = vcat( [a[:]' for a in dgrid(x, dim)] )

"compute distance of points from some point"
dist(X, X0) = sqrt(sumabs2(X .- X0, 1));
dist(X) = sqrt(sumabs2(X, 1))


"return a lattice ball"
function lattice_ball(;R=5.0, A=Atri, x0=[0.0;0.0])
    # generate a cubic portion of Zd containing -R/sig : R/sig in each
    # coordinate direction
    sig = minimum(svd(A)[2])
    ndim = ceil(Int, (norm(x0)+R)/sig)
    Z = meshgrid_list(-ndim:ndim, -ndim:ndim)
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


# simple domain constructor, just passing positions
Domain(X::Matrix{Float64}) = Domain(Matrix{Int}(), Vector{Int}(),
                     Matrix{Float64}(), size(X,2), Triangulation(X), Dict())


function Domain(; A=nothing, Ra=5.0, Rc=0.0,
                lattice=:triangular, defect=:none, shape = :ball,
                meshparams = [1.5; 3.0], x0 = :auto, V = nothing,
                edgevacancy = false, bvec = 1.0, xicorr = true )

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
      if defect == :screw || defect == :edge
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

   # shift the origin into 0
   X = X .- x0
   x0 = [0.0;0.0]

    # create a Triangulation object without actual triangulation
    # the proper triangulation will be created after we add the defects
    X = reshape(X, 2, length(X) ÷ 2)
    tri = FEM.Triangulation()
    tri.X = X

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
            add_atom!(geom, x0)
         elseif defect == :edge
            edge_predictor!(geom, AA, V; b=bvec, xicorr=xicorr)
            if edgevacancy
               Xnew = positions(geom)
               rnew = sumabs2(Xnew, 1) |> sqrt
               I0 = find(rnew .< 1)
               remove_atom!(geom, I0[2])
            end
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


"add an atom to the geometry, usually to create an interstitial"
function add_atom!(geom::Domain, x0::Vector)
   X = [positions(geom) x0]
   geom.mark = [geom.mark; 0]
   geom.nA = geom.nA + 1
   # redo the triangulation
   geom.tri = FEM.Triangulation(X)
   return geom
end



#    DISLOCATION STUFF


"compute lame parameters for given potential"
function lame_parameters(V::LujiaLt.Potentials.LennardJonesPotential)
   X, _ = lattice_ball(;R=V.cutoff[2] + 1.1)
   h = 1e-4
   dV = r -> LujiaLt.Potentials.lj1(r, V.cutoff)
   ddV = r -> (dV(r+h) - dV(r-h))/(2*h)
   μ = 0.0
   λ = 0.0
   for n = 1:size(X,2)
      r = norm(X[:,n])
      if r > 0
         μ += 0.25 * (r^2 * ddV(r) - r * dV(r))
         λ += r * dV(r)
      end
   end
   λ += μ
   ν = 0.5 * λ / (μ + λ)
   return μ, λ, ν
end


function nndist(V::LujiaLt.Potentials.LennardJonesPotential)
   X, _ = lattice_ball(;R=V.cutoff[2] + 1.1)
   rr = sumabs2(X, 1) |> sqrt
   dV = r -> r > 0 ? LujiaLt.Potentials.lj1(r, V.cutoff) : 0.0
   stress = t -> sum([ t*r * dV(t*r)  for r in rr ])
   t1 = 1.0; s1 = stress(t1)
   t2 = 0.95; s2 = stress(t2)
   while abs(s2) > 1e-8
      tnew = t2 - s2 * (t2 - t1) / (s2 - s1)
      t2, t1 = tnew, t2
      s2, s1 = stress(t2), s2
   end
   return t2
end

"isotropic CLE edge dislocation solution"
function ulin_edge_isotropic(X, b, ν)
    x, y = X[1,:], X[2,:]
    r² = x.^2 + y.^2
    r = sqrt(r²)
    ux = b/(2*π) * ( angle(x + im*y) + (x .* y) ./ (2*(1-ν) * r²) )
    uy = -b/(2*π) * ( (1-2*ν)/(4*(1-ν)) * log(r²) + - 2 * y.^2 ./ (4*(1-ν) * r²) )
    return [ux; uy]
end



"lattice corrector to CLE edge solution"
function xi_solver(Y::Vector, b; TOL = 1e-10, maxnit = 5)
   ξ1(x::Real, y::Real, b) = x - b * angle(x + im * y) / (2*π)
   dξ1(x::Real, y::Real, b) = 1 + b * y / (x^2 + y^2) / (2*π)
    y = Y[2]
    x = y
    for n = 1:maxnit
        f = ξ1(x, y, b) - Y[1]
        if abs(f) <= TOL; break; end
        x = x - f / dξ1(x, y, b)
    end
    if abs(ξ1(x, y, b) - Y[1]) > TOL
        warn("newton solver did not converge; returning input")
        return Y
    end
    return [x, y]
end

"Ehrlacher/Ortner/Shapeev edge dislocation solution"
function ulin_edge_eos(X, b, ν)
    Xmod = zeros(X)
    for n = 1:size(X,2)
        Xmod[:, n] = xi_solver(X[:,n], b)
    end
    return ulin_edge_isotropic(Xmod, b, ν)
end

# TODO: generalise this method to arbitrary potentials!
function edge_predictor!(geom::Domain, A::Matrix, V::LujiaLt.Potentials.LennardJonesPotential; b = nndist(V),
      xicorr = true)
   @assert vecnorm(A - Atri) < 1e-10
   X = positions(geom)
   X *= b
   if xicorr
      X += ulin_edge_eos(X, b, 0.25)
   else
      X += ulin_edge_isotropic(X, b, 0.25)
   end
   geom.tri = FEM.Triangulation(X)
   return geom
end
