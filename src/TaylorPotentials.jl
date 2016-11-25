
module TaylorPotentials

import LujiaLt.TightBinding: TBPotential, tb_energy1
import LujiaLt: dist
using LujiaLt.MDTools
# import LujiaLt: FrcHybridModel

export LinearSiteForce, TSiteForce, RefNeigList, Stencil

abstract TSiteForce

"""
`type Stencil`: encodes information about an interaction stencil
"""
type Stencil
   Z::Matrix{Int}
   X::Matrix{Float64}
   A::Matrix{Float64}
   rcut::Float64
   I0::Int
end


function Stencil(Xmm::Matrix{Float64}, Zmm::Matrix{Int})
   # deformation matrix: Z = A * X
   A = Xmm * pinv(Zmm)
   @assert vecnorm(A * Zmm - Xmm, Inf) < 1e-12
   # system size
   d, N = size(Xmm)
   # find the origin
   r = dist(Xmm)
   I0 = find(r .< 1e-10)[1]
   # compute nn distance and resulting finite-difference h and cut-off
   #  (we deliberate remove r[I0] from the minimum)
   r[I0] = maximum(r) + 0.001
   # rnn = minimum(r)
   # h = 1e-4 * rnn
   rcut = maximum(r)
   return Stencil(Zmm, Xmm, A, rcut, I0)
end

import Base.size
size(s::Stencil) = size(s.X)



"""
### linear site force

stencil: subset of R^d, corresponds to direction vectors in Xref
coeffs: dim x dim x (number of stencil entries)

note this actually approximates the gradient (negative force) not the
 actual force
"""
type LinearSiteForce <: TSiteForce
   stencil::Stencil
   coeffs::Array{Float64, 3}
end

LinearSiteForce(Vqm::TBPotential, Xmm::Matrix{Float64}, Zmm::Matrix{Int}) =
   LinearSiteForce(Vqm, Stencil(Xmm, Zmm))

# TODO: need to allow Yref as well here!!!!!!

LinearSiteForce(Vqm::TBPotential, stencil::Stencil) =
   LinearSiteForce(stencil, train_lin(Vqm, stencil))

cutoff(lsf::LinearSiteForce) = lsf.stencil.rcut


"""
Compute the coefficients for a linear approximation to the site force
"""
function train_lin(Vqm::TBPotential, stencil::Stencil)
   # force in reference configuration
   #    ∂_bn f_a0 = ∂_a0 f_bn
   d, N = size(stencil)
   coeffs = zeros(d, d, N)
   h = 1e-4    # TODO: this should be 1e-4 * rnn
   for a = 1:d
      stencil.X[a, stencil.I0] += h
      Fplus = tb_energy1(Vqm, stencil.X)
      stencil.X[a, stencil.I0] -= 2*h
      Fminus = tb_energy1(Vqm, stencil.X)
      coeffs[a, :, :] = (Fplus - Fminus) / (2*h)
      stencil.X[a, stencil.I0] += h
   end
   return coeffs
end



function train_neighbourlist(Xref, Vmm::TSiteForce, Imm)
   nlist = NeighbourList(Xref, cutoff(Vmm))
   Npair = length(nlist.i)
   Ncoeff = size(Vmm.stencil.Z, 2)
   icoeff = zeros(Int, Npair)
   isMM = fill(false, size(Xref, 2))
   isMM[Imm] = true
   for t = 1:Npair
      if !isMM[nlist.i[t]]
         continue
      end
      # reverse distance
      Rz = round(Int, Vmm.stencil.A \ nlist.R[:, t])
      # find the corresponding index
      # TODO: this is incredibly inefficient, it can be done better!
      for s = 1:Ncoeff
         if Rz == Vmm.stencil.Z[:, s]
            icoeff[t] = s
            break
         end
      end
   end
   Ikeep = find(icoeff .!= 0)
   return RefNeigList(nlist.i[Ikeep], nlist.j[Ikeep], icoeff[Ikeep])
end


type RefNeigList
   i::Vector{Int}
   j::Vector{Int}
   icoeff::Vector{Int}
end


# fast function to evaluate a SiteLinearForce
"""
`function mm_forces`

compute MM force contribution from `TaylorSiteForce`
"""
function mm_forces( Vmm::LinearSiteForce,
                     U::Array{Float64, 2},
                     nlist::RefNeigList )
   # extract dimension information (d1, d2 \in \{2, 3\})
   d1, d2, nneig = size(Vmm.coeffs)
   @assert size(U, 1) == d1
   # allocate output
   Frc = zeros(size(U))
   # do the actual assembly (this returns Frc again)
   return mm_frc_lin!(Vmm, U, nlist, Frc)
end


# TODO: we can probably speed this up by doing the
#        reinterpret trick
function mm_frc_lin!{T <: TSiteForce}(Vmm::T,
                     U::Array{Float64, 2},
                     nlist::RefNeigList, Frc::Matrix{Float64} )
   # loop over neighbour - pairs
   d1, d2, _= size(Vmm.coeffs)
   for (i, j, icoeff) in zip(nlist.i, nlist.j, nlist.icoeff)
      for a = 1:d1, b = 1:d2
         Frc[a, i] += Vmm.coeffs[a,b,icoeff] * (U[b,j]-U[b,i])
      end
   end
   return Frc
end


"""
`type QuadraticSiteForce`
"""
type QuadraticSiteForce <: TSiteForce
   stencil::Stencil
   coeffs::Array{Float64, 3}    # dim x dim x N
   coeffs2::Array{Float64, 5}   # dim x (dim x N) x (dim x N)
end

QuadraticSiteForce(Vqm::TBPotential, stencil::Stencil) =
   QuadraticSiteForce(stencil,
                     train_lin(Vqm, stencil),
                     train_quad(Vqm, stencil) )


"""
Compute the coefficients for a linear approximation to the site force
"""
function train_quad(Vqm::TBPotential, stencil::Stencil)
   # force in reference configuration
   #    ∂_bn f_a0 = ∂_a0 f_bn
   d, N = size(stencil)
   coeffs2 = zeros(d, d, N, d, N)
   h = 1e-4    # TODO: this should be 1e-4 * rnn
   for a=1:d, b=1:d, n=1:N
      stencil.X[a, stencil.I0] += h
      Fplus = tb_energy1(Vqm, stencil.X)
      stencil.X[a, stencil.I0] -= 2*h
      Fminus = tb_energy1(Vqm, stencil.X)
      coeffs[a, :, :] = (Fplus - Fminus) / (2*h)
      stencil.X[a, stencil.I0] += h
   end
   return coeffs
end



# # fast function to evaluate a SiteLinearForce
# function mm_forces( Vmm::QuadraticSiteForce,
#                      U::Array{Float64, 2},
#                      nlist::RefNeigList )
#    # extract dimension information (d1, d2 \in \{2, 3\})
#    d1, d2, nneig = size(Vmm.coeffs)
#    @assert size(U, 1) == d1
#    # allocate output
#    Frc = zeros(size(U))
#    # do the actual assembly (this returns Frc again)
#    return mm_frc_lin!(Vmm, U, nlist, Frc)
# end


end
