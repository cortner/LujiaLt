
import .Potentials: SitePotential, StandardSitePotential, rdim


"return a string label that described the model"
label{T <: Model}(::T) = error("no label defined for $T")
export label


################### Some generic dof-conversion functions
# any model that wants to use these must have the fields
# Ifree, Yref, V, geom
# NOTE: at the moment, these function all assume that the
#       dofs are the positions (or displacements) of the free atoms
#       if this turns out to break later models, then one can easily
#       redirect this, based on the type of the potential using, e.g.,
#   dofs2defm(m::Model, dofs::Vector{Float64}) =
#          dofs2defm(m, m.V, dofs)
#   and then dispatch on V.

ndofs(m::Model) = rdim(m.V) * length(m.Ifree)

dof_vector(m::Model) = zeros(ndofs(m))

"convert a dof-vector to a (generalised) deformation matrix"
function dofs2defm(m::Model, dofs::Vector{Float64})
    Y = copy(m.Yref)
    Y[:, m.Ifree] = reshape(dofs, rdim(m.V), length(m.Ifree))
    return Y
end

defm2dofs(m::Model, Y::Matrix{Float64}) =
    Y[:, m.Ifree][:]

"convert a force (or displacement) to a force acting on the dof-space"
frc2dofs(m::Model, P::Matrix{Float64}) =
    P[:, m.Ifree][:]

"convert a dof-type vector into a force array"
function dofs2frc(m::Model, dof::Vector{Float64})
    P = zeros(rdim(m.V), nX(m.geom))
    P[:, m.Ifree] = reshape(dof, rdim(m.V), length(m.Ifree))
    return P
end

"""returns a list of indices Jfree, where,
if y = Y[:], then y[Jfree] = Y[:, Ifree][:]
"""
function free_defm_indices(m::Model)
   J = reshape(collect(1:rdim(m.V) * nX(m.geom)),
               rdim(m.V), nX(m.geom))
   return J[:, m.Ifree][:]
end


"""
reference_configuration(geom, V)

check that geom, V are compatible and if so, return a copy of the positions
stored in geom to be stored as a reference configuration
"""
function reference_configuration(geom::Domain, V::StandardSitePotential)
    if rdim(V) != size(positions(geom), 1)
        error("reference_configuration : need rdim(V) == dDim(geom)!")
    end
    return copy(positions(geom))
end



##################### Basic Atomistic Model ###################################

"""
Atomistic (classical potential) cluster model

## Canonical Constructor:

`Atm(V::SitePotential, Ra; lattice=:triangular, defect=:vacancy)`

* Ra : core region radius of free atoms (buffer is added automatically)
* See `?Domain` for `lattice` and `defect` parameters.
"""
type Atm{TV <: SitePotential} <: Model
   geom::Domain
   V::TV
   Ifree::Vector{Int}
   Yref::Matrix{Float64}
   # ------ Atm specific fields: ----------
   vol::Vector{Float64}
end
export Atm

# default constructor for Atm
function Atm(; V=nothing, Ra=5.0, lattice=:triangular, defect=:none)
   geom = Domain(;Ra=Ra+2*cutoff(V)+0.1, lattice=lattice, defect=defect)
   r = dist(positions(geom))
   Ifree = find(r .< Ra)
   vol = zeros(nX(geom))
   vol[find(r .< Ra+cutoff(V)+0.1)] = 1.0
   Yref = reference_configuration(geom, V)
   return Atm(geom, V, Ifree, Yref, vol)
end

label(::Atm) = "ATM"

# evaluate the energy of the model
evaluate(m::Atm, dofs::Vector{Float64}) =
    at_energy_diff(m.V, dofs2defm(m, dofs), m.Yref, m.vol)

# evaluate the gradient of the energy of this model
grad(m::Atm, dofs::Vector{Float64}) =
    frc2dofs(m, at_energy1(m.V, dofs2defm(m, dofs), m.vol))


##################### B-QCE Model ###################################




##################### B-QCF Model ###################################
