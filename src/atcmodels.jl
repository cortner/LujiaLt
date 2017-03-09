

import .Potentials: SitePotential, StandardSitePotential, rdim, cutoff, nndist
import .FEM: Triangulation, elements

export BQCE, Atm, label


"return a string label that described the model"
label{T <: Model}(::T) = error("no label defined for $T")


"""
reference_configuration(geom, V)

check that geom, V are compatible and if so, return a copy of the positions
stored in geom to be stored as a reference configuration
"""
function reference_configuration(geom::Domain, V::SitePotential)
    if rdim(V) != size(positions(geom), 1)
        error("reference_configuration : need rdim(V) == dDim(geom)!")
    end
    return copy(positions(geom))
end


nndist(m::Model) = nndist(m.V)
rdim(m::Model) = rdim(m.V)
nX(m::ACModel) = nX(m.geom)

##################### Basic Atomistic Model ###################################

"""
Atomistic (classical potential) cluster model

## Canonical Constructor:

`Atm(V::SitePotential, Ra; lattice=:triangular, defect=:vacancy)`

* Ra : core region radius of free atoms (buffer is added automatically)
* See `?Domain` for `lattice` and `defect` parameters.
"""
type Atm{TV <: SitePotential} <: ACModel
   geom::Domain
   V::TV
   Ifree::Vector{Int}
   Yref::Matrix{Float64}
   fext::Matrix{Float64}
   # ------ Atm specific fields: ----------
   vol::Vector{Float64}
end

# default constructor for Atm
function Atm(; V=nothing, Ra=5.1, defect=:none, lattice=:triangular)
   geom = Domain(Ra = Ra + 2 * tight_buffer(V), defect=defect, lattice=lattice)
   r = dist(positions(geom))
   Ifree = find(r .< Ra)
   vol = zeros(nX(geom))
   vol[find(r .< Ra + tight_buffer(V))] = 1.0
   Yref = reference_configuration(geom, V)
   fext = zeros(size(Yref))
   return Atm(geom, V, Ifree, Yref, fext, vol)
end

# alternative constructor passing X directly, this creates
# a geom object without all the detailed information that the
# "proper" construction provides
Atm(X::Matrix{Float64}, Ifree::Vector{Int}, V::SitePotential) =
   Atm(Domain(X), V, Ifree, X, zeros(X), ones(size(X,2)))


label(::Atm) = "ATM"

# evaluate the energy of the model
evaluate(m::Atm, dofs::Vector{Float64}) =
    ( at_energy_diff(m.V, dofs2defm(m, dofs), m.Yref, m.vol)
      - vecdot(m.fext, dofs2defm(m, dofs) - m.Yref)   )

# evaluate the gradient of the energy of this model
grad(m::Atm, dofs::Vector{Float64}) =
    ( frc2dofs(m, at_energy1(m.V, dofs2defm(m, dofs), m.vol))
      - frc2dofs(m, m.fext) )


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
