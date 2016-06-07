

import .Potentials: SitePotential, StandardSitePotential, rdim, cutoff, nndist
import .FEM: Triangulation, elements

export BQCE, Atm, label


"return a string label that described the model"
label{T <: Model}(::T) = error("no label defined for $T")

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

Aref(m::ACModel) = m.geom.A

tight_buffer(V::SitePotential) = cutoff(V) + 0.1
moderate_buffer(V::SitePotential) = 1.1*(0.3+cutoff(V))


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

# default constructor for Atm
function Atm(; V=nothing, Ra=5.1, kwargs...)
   geom = Domain(Ra = Ra + 2 * tight_buffer(V), kwargs...)
   r = dist(positions(geom))
   Ifree = find(r .< Ra)
   vol = zeros(nX(geom))
   vol[find(r .< Ra + tight_buffer(V))] = 1.0
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

type BQCE{TV}  <: ACModel
   geom::Domain
   V::TV
   Ifree::Vector{Int}
   Yref::Matrix{Float64}
   # ------ BQCE specific fields: ----------
   volX::Vector{Float64}   # volumes associated with atoms
   volT::Vector{Float64}   # volumes associated with elements
end


"natural cubic spline on [0, 1]"
cubicspline(t) = (t .>= 1.0) + (0.0 .< t .< 1.0) .* (3*t.^2 -2*t.^3)

"standard radial blending function (cubic spline)"
radial_blending(X, r0, r1; x0=zeros(size(X,1))) =
   cubicspline( (dist(X .- x0) - r0) / (r1-r0) )[:]


"""
`BQCE(kwargs...)`

initialise a BQCE model.

* `V` : interatomic potential
* `Ra` : inner atomistic radius
* `Rb` : blending radius; note the blending width is approx. `Rb-Ra`
* `Rc` : outer (continuum) domain radius
* `lattice` : lattice type; see `?Domain`
* `defect` : lattice defect; see `?Domain`
"""
function BQCE(; V=nothing, Ra=3.1, Rb=6.1, Rc=12.1, blending=:radial,
                kwargs... )
   # compute a "good" radius for the atomistic part of the Domain
   Rab_buf = Rb + moderate_buffer(V)
   @assert Rc > Rab_buf
   # create the geometry
   geom = Domain(Ra=Rab_buf, Rc=Rc, kwargs...)
   Yref = reference_configuration(geom, V)
   # free nodes
   r = dist(positions(geom))
   Ifree = find(r .< Rc-0.1)
   # compute the blending
   # TODO: implement a displatch mechanism for blending
   @assert blending == :radial
   volX = 1.0-radial_blending(positions(geom), Ra, Rb)
   # volX defines a blending function from which we can now
   # compute the quadratue weights for the blended CB model
   # NOTE: if we want to extend to P2-fem, then we need to
   #       change this implementation quite a bit.
   volT = Float64[ el.vol / det(Atri) * mean(1.0-volX[el.t])
                   for el in elements(geom.tri) ]
   # ready to complete the assembly
   return BQCE(geom, V, Ifree, Yref, volX, volT)
end

label(::BQCE) = "BQCE"


# evaluate the energy of the model
evaluate(m::BQCE, dofs::Vector{Float64}) = evaluate(m, dofs2defm(m, dofs))
evaluate(m::BQCE, Y::Matrix{Float64}) =
               at_energy_diff(m.V, Y, m.Yref, m.volX) + cb_energy_diff(m, Y)

# evaluate the gradient of the energy of this model
grad(m::BQCE, Y::Matrix{Float64}) =
      at_energy1(m.V, Y, m.volX) + cb_energy1(m, Y)
grad(m::BQCE, dofs::Vector{Float64}) = frc2dofs(grad(m, dofs2defm(m, dofs)))


##################### B-QCF Model ###################################
