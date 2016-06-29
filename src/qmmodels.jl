
using LujiaLt.TightBinding
TB = LujiaLt.TightBinding

import LujiaLt: evaluate, grad
export QModel, TBModel

abstract QModel <: Model

"""
`type TBModel`

Encapsulates a TB model with clamped (Dirichlet) boundary conditions
"""
type TBModel <: QModel
   # geom::Domain
   Xref::Matrix{Float64}
   V::TBPotential
   Ifree::Vector{Int}
   ######
   Yref::Matrix{Float64}
   Eref::Vector{Float64}
end

function TBModel(; V=nothing, Rq=5.1, Rbuf=1+log(Rq),
                     defect=:none, lattice=:triangular)
   geom = Domain(Ra = Rq + Rbuf, defect=defect, lattice=lattice)
   Xref = positions(geom)
   r = dist(Xref)
   Ifree = find(r .<= Rq)
   Yref = TB.reference_configuration(Xref, V)  # TODO: allow `defect` argument?
   Eref = tb_energies(V, Yref)
   return TBModel(Xref, V, Ifree, Yref, Eref)
end


# evaluate the energy of the model
# we use precomputed site energies
evaluate(m::TBModel, dofs::Vector{Float64}) =
   sum_kbn( tb_energies( m.V, dofs2defm(m, dofs) ) - m.Eref )

# evaluate the gradient of the energy of this model
grad(m::TBModel, dofs::Vector{Float64}) =
    frc2dofs(m, tb_energy1(m.V, dofs2defm(m, dofs)))


nX(m::QModel) = size(m.Xref, 2)


##################### HYBRID MODELS ---- TODO ##########################


abstract HybridModel

type EnHybridModel <: HybridModel
   geom::Domain
   V::TBPotential
   Iqm::Vector{Int}
   Imm::Vector{Int}
   Ifree::Vector{Int}
   Yref::Matrix{Float64}
end

type FrcHybridModel <: HybridModel
   geom::Domain
   V::TBPotential
   Iqm::Vector{Int}
   Imm::Vector{Int}
   Ifree::Vector{Int}
   Yref::Matrix{Float64}
end
