
"""
# `module FEM`

Module that collects some basic finite element functionality for
application with the Cauchy--Born model. The core is a Julia
wrapper for Qhull which goes via `PyCall` and `scipy.spatial.Delaunay`.

## Quick example

```
using LujiaLt.FEM
tri = Triangulation(rand(2, 10))
plot(tri; xradius=0.02)
```

## List of types and methods

Look at the respective docs for more info.

* `Triangulation` : Delauany triangulation from arbitrary set of points
* `nT, nX` : number of elements / number of nodes
* `elements` : iterator over elements
* `P1_element` : some useful generic computations associated with P1-FEM
* `plot` : display the mesh
* `locate` : find elements to which some points belong

## Dependencies:

* `Compose` for plotting
* `PyCall` for calling Python / scipy
* working Python installation with `scipy`

## TODO

* Constrained Delaunay triangulation
* maybe admitting different triangulation packages
"""
module FEM

using PyCall
using Compose

@pyimport scipy.spatial as scipy_spatial

export Triangulation, locate, nT, nX
export elements, P1_element, ∇u
export plot


"""
`type Triangulation`

An elementary triangulation wrapper type, to allow in principle usage of
different background triangulation packages. The current implementation uses
`scipy.spatial.Delaunay`, which is itself a [Qhull](http://www.qhull.org)
wrapper.

* `X` : dim x Npoints array containing the vertices
* `T` : dim+1 x Ntriangles array containing the element information

## Constructor

* To construct a Delaunay triangulation create a 2 x Npoints or a Npoints x 2
array `X`, where Npoints must be larger than 2, and simply call
`tri = Triangulation(X)`.

## Usage via iterator

TODO: update documentation here!

```
  # assume X are the FE nodes and U the nodal values
  tri = Triangulation(X)
  for iT in elements(tri)
      vol, F, Du = P1_element(iT, tri, U)
      # do something with it
  end
```
"""
type Triangulation
    X::Matrix{Float64}
    T::Matrix{Int}
    pyo::PyObject
end

"""
`type P1element`

returned during iterations with `elements(tri)`.

## Fields

* `t`: indices of vertices
* `vol` : volume of element
* `J` : jacobian matrix for coordinate transformation
* `B` : P1 gradient operator: B * U = ∇u (see also `?∇u`)
* `xT` : element mid-point
* `idx` : element index; but note that `A[el.idx]` is equivalent to `A[el]`
"""
type P1element
  t::Vector{Int}
  vol::Float64
  J::Matrix{Float64}
  B::Matrix{Float64}
  xT::Vector{Float64}
  idx::Int
end

# a cute getindex hack
Base.getindex(coll::AbstractArray, el::P1element) = coll[el.idx]
Base.setindex!(coll::AbstractArray, x, el::P1element) = (coll[el.idx] = x)

# empty triangulation constructor
Triangulation() = Triangulation(Matrix{Float64}(), Matrix{Int}(), PyObject(nothing))

function Triangulation(X)
    # reshape X
    if size(X, 1) == 2
        X = X'
    elseif size(X, 2) != 2
        error("Triangulation : size(X) must be (2, N) or (N, 2)")
    end
    # call the Qhull code via PyCall via scipy bindings
    tri = scipy_spatial.Delaunay(X)
    # return the Julia object
    return Triangulation(X', tri[:simplices]' + 1, tri)
end

"return number of elements in the triangulation"
nT(tri::Triangulation) = size(tri.T, 2)

"returns number of nodes in the triangulation"
nX(tri::Triangulation) = size(tri.X, 2)

"creates an iterator over elements"
elements(tri::Triangulation) = tri
Base.start(::Triangulation) = 0
Base.done(tri::Triangulation, state::Int) = (state == nT(tri))
Base.next(tri::Triangulation, state::Int) = ( P1element(state+1, tri), state+1 )
Base.length(tri::Triangulation) = nT(tri)

# TODO: consider reviving this!
# TriU = Tuple{Triangulation, Vector}
# elements(tri::Triangulation, U) = (tri, U)
# Base.start(::TriU) = 0
# Base.done(triU::TriU, state::Int) = (state == nT(triU[1]))
# Base.next(triU::TriU, state::Int) = ( P1element(state+1, triU...), state+1 )


"compute the P1 gradient in the reference triangle"
ref_grad(t::Vector{Int}, V::Matrix{Float64}) = V[:, t[2:3]] .- V[:, t[1]]

∇u(el, U::Vector) = el.B' * U[el.t]
∇u(el, U::Matrix) = U[:, el.t] * el.B

stress2frc(el, T::Matrix) = T * el.B'

"""
compute some information related to P1-FEM:

* `P1element(idx::Int, tri::Triangulation) -> el::P1element`
* `P1_element(el::Int, tri::Triangulation, U) -> vol, F, Du, u, xT`,
    where `Du` of dimension r x dim (dim=2) is the gradient of the P1 function
        represented by `U` which is given as a `r x nX` array with `r` being the
        range dimension and `u` is the value of `u` at the element mid-point.
"""
function P1element(idx::Int, tri::Triangulation)
   t = tri.T[:, idx]
   F = ref_grad(t, tri.X)
   return P1element( tri.T[:, idx], det(F) / 2, F,
                      [ -1 -1; 1 0; 0 1] / F, mean(tri.X[:, t], 2)[:], idx )
end

# function P1element(idx::Int, tri::Triangulation, U::Matrix{Float64})
#    el = P1element(idx, tri)
#    return el, el.B' * U[el.t]
# end
# P1element( eltri::Tuple, varargs...) = P1_element(eltri..., varargs...)

"""
return the convex coordinates of a point `y` relative to a triangle
with index `i`
"""
function convex_coordinates(y::Vector, i::Integer, tri::Triangulation)
   # y = x0 + F * [l1;l2]
   F = ref_grad(tri.T[:, i], tri.X)
   lam = F \ (y - tri.X[:, tri.T[1,i]])
   return [1 - lam[1] - lam[2]; lam[1]; lam[2]]
end


"""
find element indices to which the provided points belong.

## Parameters
* `P` : 2 x M, list of points that we want to locate
* `tri` : Triangulation
"""
locate(P::Matrix{Float64}, tri::Triangulation) = tri.pyo[:find_simplex](P') + 1


"""
compute unique set of edges: for now this is an auxiliary function used
for plotting, but it could become something more useful for other
purposes as well
"""
function edges(tri::Triangulation)
    T = tri.T
    S = [ T[1,:] T[1,:] T[2,:]; T[2,:] T[3,:] T[3,:] ]
    S = sort(S, 1)
    return unique(S, 2)
end


end # module FEM
