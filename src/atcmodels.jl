
import .Potentials: SitePotential, StandardSitePotential, rdim

export Atm, label

"return a string label that described the model"
label{T <: Model}(::T) = error("no label defined for $T")


##################### Basic Atomistic Model ###################################

"""
Atomistic (classical potential) cluster model

## Canonical Constructor:

`Atm(V::SitePotential, Ra; lattice=:triangular, defect=:vacancy)`

* Ra : core region radius of free atoms (buffer is added automatically)
* See `?Domain` for `lattice` and `defect` parameters.
"""
type Atm{T <: SitePotential} <: Model
    geom::Domain
    V::T
    Ifree::Vector{Int}
    vol::Vector{Float64}
    Yref::Matrix{Float64}
end

function reference_configuration(geom::Domain, V::StandardSitePotential)
    if rdim(V) != size(positions(geom), 1)
        error("reference_configuration : need rdim(V) == dDim(geom)!")
    end
    return copy(positions(geom))
end

function Atm(; V=nothing, Ra=0.0, lattice=:triangular, defect=:none)
    geom = Domain(;Ra=Ra, Rbuf=2*cutoff(V)+0.1, lattice=lattice, defect=defect)
    r = dist(positions(geom))
    Ifree = find(r .< Ra)
    vol = zeros(nX(geom))
    vol[find(r .< Ra+cutoff(V)+0.1)] = 1.0
    Yref = reference_configuration(geom, V)
    return Atm(geom, V, Ifree, vol, Yref)
end

label(::Atm) = "ATM"
    
dof_vector{T <: StandardSitePotential}(m::Atm{T}) = zeros(rdim(m.V) * length(m.Ifree))

"convert a dof-vector to a (generalised) deformation matrix"
function dofs2defm{T <: StandardSitePotential}(m::Atm{T}, dofs::Vector{Float64})
    Y = copy(m.Yref)
    Y[:, m.Ifree] = reshape(dofs, rdim(m.V), length(m.Ifree))
    return Y
end

defm2dofs{T <: StandardSitePotential}(m::Atm{T}, Y::Matrix{Float64}) =
    Y[:, m.Ifree][:]

"convert a force (or displacement) to a force acting on the dof-space"
frc2dofs{T <: StandardSitePotential}(m::Atm{T}, P::Matrix{Float64}) = 
    P[:, m.Ifree][:]

"convert a dof-type vector into a force array"
function dofs2frc{T <: StandardSitePotential}(m::Atm{T}, dof::Vector{Float64})
    P = zeros(rdim(m.V), nX(m.geom))
    P[:, m.Ifree] = reshape(dof, rdim(m.V), length(m.Ifree))
    return P
end

# evaluate the energy of the model
evaluate(m::Atm, dofs::Vector{Float64}) =
    at_energy(m.V, dofs2defm(m, dofs), m.vol)

# evaluate the gradient of the energy of this model
grad(m::Atm, dofs::Vector{Float64}) =
    frc2dofs(m, at_energy1(m.V, dofs2defm(m, dofs), m.vol))


##################### B-QCE Model ###################################




##################### B-QCF Model ###################################

