
"""
# `module FEM`

Module that collects some basic finite element functionality for
application with the Cauchy--Born model. 

## List of types

Look at the respective documentions for more info.

* `Triangulation` : Delauany triangulation from arbitrary set of points
* `nT, nX` : number of elements / number of nodes
* `elements` : iterator over elements
* `P1_element` : some useful generic computations associated with P1-FEM
* `plot` : display the mesh

## Dependencies:

* Gadfly for plotting 
* PyCall for calling Python Triangle
* [Python Triangle](https://github.com/drufat/triangle)

"""
module FEM

using PyCall
using Compose

@pyimport triangle

export Triangulation
export nT, nX
export elements
export P1_element
export plot

"""`type Triangulation`

An elementary triangulation wrapper type, to allow in principle usage of
different background triangulation packages. The current implementation
requires [Python Triangle](https://github.com/drufat/triangle)

* `X` : dim x Npoints array containing the vertices
* `T` : dim+1 x Ntriangles array containing the element information

## Constructor

* To construct a Delaunay triangulation create a 2 x Npoints or a Npoints x 2
array `X`, where Npoints must be larger than 2, and simply call
`tri = Triangulation(X)`.

* For now the code does not allow constrained Delaunay. If the domain is
just barely convex (e.g. a polygon), then this means that thin elements
will be created at the boundary. To account for that, if a list of 
boundary node indices is passed as well, then the domain will first be
artifically convexified before triangulation:
```
  tri = Triangulation(X; bdry=list_of_indices)
```

## Usage via iterator

```
# assume X are the FE nodes an U the nodal values
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
end


function Triangulation(X; bdry=nothing)
    # convexify if needed
    if bdry != nothing
        # do the convexification!
    else
        Y = X
    end
    # reshape X
    if size(Y, 1) == 2
        Y = Y'
    elseif size(Y, 2) != 2
        error("Triangulation : size(X) must be (2, N) or (N, 2)")
    end
    # call the Shewchuck-code via python
    tri = triangle.triangulate(Dict("vertices" => Y))
    # return the Julia object
    return Triangulation(Y', tri["triangles"]' + 1)
end

"return number of elements in the triangulation"
nT(tri::Triangulation) = size(tri.T, 2)

"returns number of nodes in the triangulation"
nX(tri::Triangulation) = size(tri.X, 2)

"auxiliary data-type to iterate over finite elements"
type ElementIterator
    tri::Triangulation
    i::Int
end

"creates an iterator over elements"
elements(tri::Triangulation) = ElementIterator(tri, 0)
done(tri::Triangulation, state::ElementIterator) = (state.i == nT(tri))
next(tri::Triangulation, state::ElementIterator) =
    (state.i+1, ElementIterator(tri, state.i+1))


ref_grad(t::Vector{Int}, V::Matrix{Float64}) =
    V[:, t[2:3]] .- V[:, t[1]]

"""compute some information related to P1-FEM:

* `P1_element(el::Int, tri::Triangulation) -> vol, F`, where `vol` is the 
    element volume and `F` the deformation matrix.
* `P1_element(el::Int, tri::Triangulation, U) -> vol, F, Du`, 
    where `Du` of dimension r x dim (dim=2) is the gradient of the P1 function
        represented by `U` which is given as a `r x nX` array with `r` being the 
        range dimension.
"""
function P1_element(el::Int, tri::Triangulation)
    F = ref_grad(tri.T[:, el], tri.X)
    return det(F) / 2, F
end

function P1_element(el::Int, tri::Triangulation, U::Matrix{Float64})
    vol, F = P1_element(el, tri)
    DU_ref = ref_grad(tri.T[:, el], U)
    invF = [F[2,2] -F[1,2]; -F[2,1] F[1,1]] * 2 / vol
    return vol, F, DU_ref * invF
end




############### PLOTTING ###############

"""compute unique set of edges: for now this is an auxiliary function use
for plotting, but it could become something more useful for 
"""
function _edges_(tri::Triangulation)
    T = tri.T
    S = [ T[1,:] T[1,:] T[2,:] T[2,:] T[3,:] T[3,:];
          T[2,:] T[3,:] T[1,:] T[3,:] T[1,:] T[2,:] ]
    S = sort(S, 1)
    return unique(S, 2)
end

# function plot(tri::Triangulation)
#     Gadfly.plot(x=tri.X[1,:], y=tri.X[2,:])
# end

function plot(tri; width=15cm, height=:auto, xradius=0.25, atcol="tomato", 
              buffer=2*xradius, elcol = "aliceblue", linecol="darkblue",
              lwidth=1.0)
    
    # create a canvas
    xLim = [extrema(tri.X[1,:])...]
    dat_width = xLim[2]-xLim[1]
    yLim = [extrema(tri.X[2,:])...]
    dat_height = yLim[2]-yLim[1]
    xLim[1] -= buffer; xLim[2] += buffer; dat_width = xLim[2]-xLim[1]
    yLim[1] -= buffer; yLim[2] += buffer; dat_height = yLim[2]-yLim[1]

    if height == :auto
        height = width * (dat_height / dat_width)
    end
    ub = UnitBox(xLim[1], yLim[1], dat_width, dat_height)
    
    # draw the segments
#     E = FEM._edges_(tri)
#     c_ctx = context(units=ub)
#     for n = 1:size(E, 2)
#         c_ctx = compoe( c_ctx, line([ tuple(tri.X[:,E[1,n]]...), tuple(tri.X[:,E[2,n]]...) ]) )
#     end
#     c_ctx = compose(c_ctx, stroke("cadetblue") ) # , linewidth(1mm))
    
    c_ctx = context(units=ub)
    for n = 1:FEM.nT(tri)
        p = tri.X[:, tri.T[:, n]]
        c_ctx = compose(c_ctx, polygon( [tuple(p[:,1]...), tuple(p[:,2]...), tuple(p[:,3]...)] ))
    end
    c_ctx = compose(c_ctx, fill(elcol), stroke(linecol), linewidth(lwidth))
    
    # draw the atoms
    a_ctx = context(units=ub)
    a_ctx = compose(a_ctx, circle(tri.X[1,:][:], tri.X[2,:][:], xradius * ones(FEM.nX(tri))), 
                    fill(atcol), stroke(linecol) )
    
    return compose(context(), a_ctx, c_ctx)
end    

end

