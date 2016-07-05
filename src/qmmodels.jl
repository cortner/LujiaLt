
using LujiaLt.TightBinding
TB = LujiaLt.TightBinding

import LujiaLt: evaluate, grad
import LujiaLt.TaylorPotentials: TSiteForce, LinearSiteForce, RefNeigList

export TBModel, FrcHybridModel

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


grad(m::TBModel, Y::Matrix{Float64}) = tb_energy1(m.V, Y)

# evaluate the gradient of the energy of this model
grad(m::TBModel, dofs::Vector{Float64}) =
    frc2dofs(m, tb_energy1(m.V, dofs2defm(m, dofs)))


nX(m::QModel) = size(m.Xref, 2)


##################### HYBRID MODELS ---- TODO ##########################



"""
`type FrcHybridModel`
"""
type FrcHybridModel{TTB <: TBPotential, TSF <: TSiteForce} <: HybridModel
   Xref::Matrix{Float64}
   V::TTB
   Ifree::Vector{Int}
   ###### predictor configuration
   Yref::Matrix{Float64}
   ###### QM/MM stuff
   Iqm::Vector{Int}
   Imm::Vector{Int}
   Vmm::TSF
   ###### neighbourlist for reference configuration
   ###### (to assemble Vmm)
   refnlist::RefNeigList
   # i::Vector{Int}
   # j::Vector{Int}
   # icoeff::Vector{Int}
end


function FrcHybridModel(; Vqm=nothing, Rqm=5.1, Rmm=20.1, Rbuf=1.0+log(Rqm),
                        Rbufmm = Rbuf,
                        defect=:none, lattice=:triangular,
                        VMM = TaylorPotentials.LinearSiteForce)

   @assert Rmm > Rqm + Rbuf + 1.0
   geom = Domain(Ra = Rmm, defect=defect, lattice=lattice)
   Xref = positions(geom)
   Yref = TB.reference_configuration(Xref, Vqm)  # TODO: allow `defect` argument?
   r = dist(Xref)
   Ifree = find( r .< Rmm - Rbuf )
   Iqm = find(r .<= Rqm)
   Imm = setdiff(Ifree, Iqm)

   # construct a little cluster to learn from
   geom_mm = Domain(Ra = Rbufmm, defect=:none, lattice=lattice)
   Xmm = positions(geom_mm)
   Zmm = convert(Matrix{Int}, geom_mm.Z)
   Vmm = VMM(Vqm, Xmm, Zmm)

   # construct the interaction information
   refnlist = TaylorPotentials.train_neighbourlist(Xref, Vmm, Imm)

   @show typeof(refnlist)

   # construct and return the thing
   return FrcHybridModel(Xref, Vqm, Ifree, Yref, Iqm, Imm, Vmm, refnlist)
end


grad(m::FrcHybridModel, Y) =
   TaylorPotentials.mm_forces(m.Vmm, Y - m.Yref, m.refnlist)



# type EnHybridModel <: HybridModel
#    geom::Domain
#    V::TBPotential
#    Iqm::Vector{Int}
#    Imm::Vector{Int}
#    Ifree::Vector{Int}
#    Yref::Matrix{Float64}
# end
