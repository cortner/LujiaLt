


module TightBinding


using LujiaLt.Potentials

import LujiaLt.Potentials: fcut, fcut1, morse, morse1 # swcut, swcut1,
import LujiaLt: ColumnIterator, columns, evaluate, grad

export TBHamiltonian, TBToyHamiltonian, TBModel, TBToyModel

export hop, hop1


"""
Abstract supertype for 2-centre tight-binding hamiltonians with only
a hoping function, but no on-site term
"""
abstract TBHamiltonian

"""
Tight-binding in-plane toy-model. No on-site term, and hopping term is
given by the morse potential.
"""
type TBToyHamiltonian
   A::Float64
   cutoff::NTuple{2, Float64}
end

hop(r, H::TBToyHamiltonian) = morse(r, H.A, H.cutoff)
hop1(r, H::TBToyHamiltonian) = morse1(r, H.A, H.cutoff)
cutoff(H::TBToyHamiltonian) = H.cutoff[2]

type TBModel{TH <: TBHamiltonian}
   H::TH   # hamiltonian
   mu::Float64      # chemical potential
   beta::Float64    # temperature / smearing parameter
   # f::Function      # smearing function
   # f1::Function     #     ... and its derivative
end


TBToyModel(; A=3.0, cutoff=(1.5, 2.3), beta = 0.1, mu = 0.0 ) =
   TBModel( TBToyHamiltonian(A, cutoff), mu, beta )


function evaluate( H::TBHamiltonian, X::Matrix{Float64} )
   d, N = size(X)
   nlist = NeighbourList(X, cutoff(H))
   H = zeros(N, N)
   for (i, j, r) in zip(nlist.i, nlist.j, nlist.r)
      H[i,j] = hop(r, H)
   end
   return H
end


function grad( H::TBHamiltonian, X::Matrix{Float64} )
   d, N = size(X)
   nlist = NeighbourList(X, cutoff(H))
   dH = zeros(d, N, N)
   nlR = reinterpret(nlist.R, )
   for (i, j, r, R) in zip(nlist.i, nlist.j, nlist.r, columns(nlist.R))
      dH[:,i,j] = (hop1(r, H)/r) * R
   end
   return dH
end

function eval_and_grad( H::TBHamiltonian, X::Matrix{Float64} )
   d, N = size(X)
   nlist = NeighbourList(X, cutoff(H))
   H = zeros(N, N)
   dH = zeros(d, N, N)
   for (i, j, r, R) in zip(nlist.i, nlist.j, nlist.r, columns(nlist.R))
      H[i,j] = hop(r, H)
      dH[:,i,j] = (hop1(r, H)/r) * R
   end
   return H, dH
end


"""
helper function to compute eigenvalues, then sort them
"""
function sorted_eig(H)
    @assert any(isnan(H)) || any(isinf(H))
    epsn, C = eig(Symmetric(H))
    Isort = sortperm(epsn)
    return epsn[Isort], C[:, Isort]
end

"""
Fermi-Dirac distribution function
"""
fermidirac( tbm::TBModel, epsn ) = 1.0 ./ (1.0 + exp(tbm.beta * (epsn-tbm.mu)))

"""
derivative of `fermidirac`
"""
fermidirac1( tbm::TBModel, epsn ) =
      - (2.0 * tbm.beta) ./ (1.0 + exp(tbm.beta * (epsn-tbm.mu))).^2 .*
      exp(tbm.beta * (epsn-tbm.mu))

"""
tight-binding site-energies
"""
function tb_energies( tbm::TBModel, X::Matrix{Float64} )
   H = evaluate(tbm.H, X)
   epsn, C = sorted_eig(H)
   f = fermidirac( tbm, epsn )
   # E_n = \sum_s ϵ_s f_s [ψ_s]_n^2    =>    (E_n) = C.^2 * (epsn .* f)
   @simd for n = 1:length(C)
      @inbounds C[n] = C[n]*C[n]
   end
   return C * (epsn .* f)
end


"""
tight-binding energy-difference
"""
tb_energy_diff( tbm::TBModel, X, Xref, vol ) =
   kbn_sum( (tb_energies( tbm, X ) - tb_energies( tbm, Xref ) ) .* vol )

tb_energy_diff( tbm::TBModel, X, Xref ) =
   tb_energy_diff( tbm, X, Xref, ones(size(X,2)) )


"""
tight-binding gradient / neg. forces
"""
function tb_energy1(tbm::TBModel, X::Matrix{Float64})
   H, dH = eval_and_grad(tbm.H, X)
   epsn, C = sorted_eig(H)
   # derivative of ϵ_s f_s
   f = fermidirac( tbm, epsn )
   df = fermidirac1( tbm, epsn )
   dg = f + epsn .* df
   # forces
   dE = 2.0 * Float64[ dot( dg, C[i,:][:] .* (C' * dH[a,i,:][:]) )
                      for a = 1:d, i = 1:N ]
   return dE
end


"""
return the density matrix (this is a very inefficient, and possibly
   numerically unstable implementation)
"""
function density_matrix( tbm::TBModel, X::Matrix{Float64} )
   H = evaluate(tbm.H, X)
   epsn, C = sorted_eig(H)
   f = fermidirac( tbm, epsn )
   N = size(X, 2)
   rho = zeros(N,N)
   for s = 1:N
      rho += f[s] * C[:, s] * C[:, s]'
   end
   return rho
end


tb_site_energy1( tbm::TBModel, X::Matrix{Float64}, I::Int ) =
      tb_site_energy( tbm, X, [I;] )


"""
`tb_site_energy1( tbm::TBModel, X::Matrix{Float64}, I::Vector{Int} )`

This is the efficient implementation of the TB site energy derivative.
"""
function tb_site_energy1( tbm::TBModel, X::Matrix{Float64}, I::Vector{Int} )

   # preparations
   H, dH = eval_and_grad(tbm.H, X)
   epsn, C = sorted_eig(H)
   f = fermidirac( tbm, epsn )
   df = fermidirac1( tbm, epsn )
   dg = f + epsn .* df

   # allocate output   >>>>  ????????
   dEs = zeros(Float64, d, N)

   # allocate temporary arrays that enable later computation
   # via in-place BLAS
   G = zeros(Float64, N, N)
   G1 = zeros(Float64, N, N)
   G2 = zeros(Float64, N, N)
   g = zeros(Float64, N)
   dHa = zeros(Float64, N, N)
   c = zeros(Float64, N)
   diff_eps_inv = zeros(Float64, N)
   D_epsn_s = zeros(Float64, N)
   if length(I) == 1
      cI = squeeze(C[I, :], 1)
   else
      cI = C[I, :]
   end
   csI_j_1 = zeros(Float64, N)


   # outer loop over dimension index - a, and over energy-level s
   for a = 1:d
      # copy a sub-array of dH for faster access
      copy!(dHa, slice(dH, a, :, :))
      # a little precomputation - BLAS-ised for performance
      Base.LinAlg.BLAS.gemm!('N', 'N', 1.0, dHa, C, 0.0, G1)
      # gs[j] = dHa[j,:] . C[:,s]

      for s = 1:N
           # copy a part of C that we need to use with BLAS
           for iii = 1:N; c[iii] = C[iii,s]; end

           # compute the G array
           scale!(dHa, c)
           G = dHa
         #   for ig = 1:N
         #       @simd for j=1:N
         #           @inbounds G[j, ig] = dHa[j,ig] * c[j]
         #       end
         #   end
           @simd for j = 1:N
               @inbounds G[j,j] += G1[j,s]
           end

           # compute the epsn_s_j > j-th entry of D_epsn_s
           # D_epsn_s =  G * c # C[:, s]
           Base.LinAlg.BLAS.gemv!('N', 1.0, G, c, 0.0,  D_epsn_s)

           # the following line was the bottle-neck at some point,
           # after moving to BLAS, not anymore
           # g = - (C' * g) ./ (epsilon - epsilon[s])
           # G2 = - G * C
           Base.LinAlg.BLAS.gemm!('N', 'N', -1.0, G, C, 0.0, G2)

           # prepare to divide by energy
           @inbounds for iii = 1:N
               diff_eps_inv[iii] = 1.0 / (epsilon[iii] - epsilon[s])
           end
           # HACK to fix the multiple e-val issue!
           for iii = max(1,s-5):min(N, s+5)
               if abs(diff_eps_inv[iii]) > 1e10
                   diff_eps_inv[iii] = 0.0
               end
           end
           # now divide by energy
           scale!(diff_eps_inv, G2)
         #   for ig = 1:N
         #       @simd for iii = 1:N
         #           @inbounds G2[iii,ig] *= diff_eps_inv[ig]
         #       end
         #   end

           # invert coordinate transform
           if length(I) == 1
               Base.LinAlg.BLAS.gemv!('N', 1.0, G2, cI, 0.0, csI_j_1)
           else
               csI_j = C[I,:] * G2'
           end

           # add the computed values to the site energy derivative
           # >>> idealised O(1)
           #  NOTE the weird construction of [C[I, s]] is just to
           #    circumvent the strange situation that if I has just
           #    one entry, then C[I,s] is a scalar while  csI_j is
           #    an array, and `dot` does not like this
           if length(I) == 1
               for j = 1:N
                   dEs[a, j] += dFermi[s] * D_epsn_s[j] * cI[s]^2 +
                                  2. * Fermi[s] * cI[s] * csI_j_1[j]
               end
           else
               for j = 1:N
                   dEs[a, j] += dFermi[s] * D_epsn_s[j] * sum( C[I,s].^2 ) +
                                  2. * Fermi[s] * dot( [C[I, s]], csI_j[:, j] )
               end
           end
      end # for s
   end # for a

   # return computed site energies
   return dEs
end


end # module TightBinding
