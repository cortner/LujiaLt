
"""
# Lujia-Light

"""
module LujiaLt

# basic molecular dynamics stuff: in particular a neighbourlist
include("mdtools.jl")

# some simple mad-up interatomic potentials
include("potentials.jl")

# some basic FEM methodology
include("fem.jl")

# this type is not in geom.jl because Plotting depends on it
# while geom again depends on Plotting.
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
export Domain


# sub-modules
include("plotting.jl")

# parts of main-module
include("geom.jl")

# assembly routines
include("assemble.jl")

# codes to help test the package
include("testing.jl")


end
