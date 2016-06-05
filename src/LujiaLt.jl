
"""
# Lujia-Light

"""
module LujiaLt

# basic molecular dynamics stuff: in particular a neighbourlist
include("mdtools.jl")

# some simple mad-up interatomic potentials
include("potentials.jl")
import .Potentials.SitePotential
# , evaluate, grad
# export evaluate, grad

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

## Constructor `Domain(kwargs...)`

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
type Domain
    # X::Matrix{Float64}
    Z::Matrix{Int32}
    mark::Vector{Int8}
    A::Matrix{Float64}
    nA::Int
    tri::FEM.Triangulation
    info::Dict
end
export Domain


"""
A `AbstractModel` is a complete description of the computational task, including
the interatomic potential (SitePotential or TB-Hamiltonian), the geometry
(Domain) whether it is fully atomistic, a/c, qm/mm, etc. The following models
have been implemented:

* ATM
* BQCE (todo)
* BQCF (todo)
* TB (todo)
* QMMM (todo)

## Computing energy and forces

If `m::Model`, and `x` a dof-vector or position array, then the call `m(x)` will
be forwarded to `evaluate(m, x)` which will compute its energy. Similarly, the
call `@GRAD m(x)` will be forwarded to `grad(m, x)` which computes the gradient
of the energy.

If the model has no associated energy, then `m(x)` should throw an error (since
it expects the energy) while `@GRAD m(x)` should return the negative force
vector.
"""
abstract Model
export Model


# some auxliary plotting functionality (needs to be defined after Domain)
include("plotting.jl")

# parts of main-module
include("geom.jl")

# some fun macros for a nicer notation
include("callmagic.jl")

# assembly routines
include("assemble.jl")

# concrete implementations of `Model` of a/c type
include("atcmodels.jl")

# interface to nonlinear solvers
include("solve.jl")

# codes to help test the package
include("testing.jl")


end
