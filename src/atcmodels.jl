
import .Potentials: SitePotential, StandardSitePotential, rdim, cutoff
using FEM

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
   J = reshape(collect(1:rdim(m.V) * nX(m.geom)), rdim(m.V), nX(m.geom))
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


positions(m) = positions(m.geom)


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



##################### General Cauchy--Born Material ###########################

abstract ACModel <: Model

"""
`getR_cb(Aref::Matrix, V::SitePotential)`

get an interaction stencil to compute Wcb with; it takes a crude estimate
on what is needed, using the assumptions (1) that interaction decays
exponentially and (2) that for large atomistic region (i.e. high accuracy),
the gradients will be close to reference.

## Parameters

* `Aref` : lattice matrix
* `V` : interaction potential
"""
function getRcb(Aref::Matrix, V::SitePotential)
   X, _ = lattice_ball(R = (nndist(V) + 0.7)*1.1, A=Aref)
   return X[:, find(sumabs2(X, 1) .> 1e-10)]
end


"Cauchy-Born energy density"
Wcb(F, V::SitePotential, R) = evaluate(V, F*R)

"derivative of Cauchy--Born energy density"
function Wcb1(F, V::SitePotential, R)
   dV = grad(V, F*R)
   dW = zeros(size(F))
   for n = 1:size(R,2)
      dW += dV[:, n] * R[:, n]'
   end
   return dW
end

"Cauchy-Born potential energy, as array of local contributions"
function cb_energies(m::ACModel, Y)
   Rcb = getRcb(m.Aref, m.V)
   Ec = zeros(nT(m.geom.tri))
   for el in elements(tri)
      if vols[el.idx] > 0
         F = ∇u(el, Y)
         Ec[el.idx] = Wcb(F, m.V, Rcb)
      end
   end
   return Ec
end

"compute energy difference between two states"
cb_energy_diff(m::ACModel, Y) = cb_energy_diff(m::ACModel, Y, m.Yref)
cb_energy_diff(m:ACModel, Y, Yref) =
                  sum_kbn( cb_energies(m, Y) - cb_energies(m, Yref) )

"Cauchy--Born energy gradient"
function cb_grad(m::ACModel, Y)
   Rcb = getRcb(m.Aref, m.V)
   dE = zeros(size(Y))
   for el in elements(tri)
      if vols[el.idx] > 0
         dE[el.t] += vol[el.idx] * el.B * Wcb1(∇u(el, Y), m.V, Rcb)
      end
   end
   return dE
end


##################### B-QCE Model ###################################




##################### B-QCF Model ###################################
