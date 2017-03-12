
import .Potentials: SitePotential, cutoff, evaluate, grad
import .MDTools: NeighbourList, sites
import .FEM: nT, ∇u, stress2frc

# export Model, Atm

########################## some generic stuff #################################

positions(m) = positions(m.geom)

Aref(m::ACModel) = m.geom.A

tight_buffer(V::SitePotential) = cutoff(V) + 0.1
moderate_buffer(V::SitePotential) = 1.1*(0.3+cutoff(V))
conservative_buffer(V::SitePotential) = (cutoff(V) + 0.7)*1.1


########################## atomistic assembly #################################

"""
* `V::SitePotential` : site energy encoding the atomistic model
* `X::Matrix{Float64}` : positions
"""
function at_energies(V::SitePotential, X::Matrix, vol)
    nlist = NeighbourList(X, cutoff(V))
    E = zeros(size(X, 2))
    for (n, neigs, r, R) in sites(nlist)
        if vol[n] > 0.0    # skip vol-0 sites
           E[n] = evaluate(V, r, R)
        end
    end
    return E
end

"""
`at_energy(V::SitePotential, X::Matrix, vol::Vector)`

compute potential energy ∑ vol[i] V(Dx(i))
"""
at_energy(V::SitePotential, X::Matrix, vol::Vector) =
   sum_kbn( vol .* at_energies(V, X, vol) )

"""
`at_energy_diff(V::SitePotential, X::Matrix, Xref::Matrix, vol::Vector)`

compute potential energy difference ∑ vol (V(Dx) - V(Dx_ref))
"""
at_energy_diff(V::SitePotential, X::Matrix, Xref::Matrix, vol::Vector) =
   sum_kbn( vol .* (at_energies(V, X, vol) - at_energies(V, Xref, vol)) )


"""
* `V::SitePotential` : site energy encoding the atomistic model
* `X::Matrix{Float64}` : positions
* `vol` : effective volumes (default 1)
"""
function at_energy1(V::SitePotential, X::Matrix{Float64}, vol::Vector{Float64})
    nlist = NeighbourList(X, cutoff(V))
    dE = zeros(X)
    for (n, neigs, r, R) in sites(nlist)
        if vol[n] <= 0.0; continue; end
        dV = vol[n] * grad(V, r, R)
        dE[:, neigs] += dV
        dE[:, n] -= sum(dV, 2)
    end
    return dE
end

# helpers to call energy/forces without passing volumes
at_energy(V::SitePotential, X::Matrix{Float64}) =
    at_energy(V, X, ones(size(X, 2)))
at_energy1(V::SitePotential, X::Matrix{Float64}) =
    at_energy1(V, X, ones(size(X, 2)))


# using LujiaLt.Potentials: LennardJonesPotential, lj1
#
# # a faster assembly for the Lennard-Jones potential
# function LujiaLt.at_energy1(V::LennardJonesPotential, X::Matrix{Float64}, vol::Vector{Float64})
#  nlist = NeighbourList(X, cutoff(V))
#  dE = zeros(X)
#  rdim = size(X,1)
#  for (n, neigs, r, R) in sites(nlist)
#     if vol[n] <= 0.0; continue; end
#     for (j, r, R) in zip(neigs, r, R)
#        dV = vol[n] * lj1(r, V.cutoff)
#        for a = 1:rdim
#           dE[a, j] += dV
#           dE[a, n] -= dV
#        end
#     end
#  end
#  return dE
# end



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
 X, _ = lattice_ball(R = conservative_buffer(V), A=Aref)
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
function cb_energies(m::ACModel, Y::Matrix, vols::Vector)
   Rcb = getRcb(Aref(m), m.V)
   Ec = zeros(nT(m.geom.tri))
   for el in elements(m.geom.tri)
      # if el.idx == 1
      #    # @code_warntype ∇u(el, Y)
      #    F = ∇u(el, Y)
      #    FR = F * Rcb
      #    r = sqrt(sumabs2(FR, 1))[:]
      #    # @code_warntype Wcb(F, m.V, Rcb)
      #    @code_warntype evaluate(m.V, r, FR)
      #    error("stop")
      # end
      if vols[el] > 0
         Ec[el] = Wcb(∇u(el, Y), m.V, Rcb)
      end
   end
   return Ec
end

cb_energies(m::ACModel, Y::Matrix) = cb_energies(m, Y, m.volT)

"compute energy difference between two states"
cb_energy_diff(m::ACModel, Y) = cb_energy_diff(m::ACModel, Y, m.Yref)
cb_energy_diff(m::ACModel, Y, Yref) =
                sum_kbn( m.volT .* (cb_energies(m, Y) - cb_energies(m, Yref)) )

"Cauchy--Born energy gradient"
function cb_energy1(m::ACModel, Y::Matrix, vols::Vector)
   Rcb = getRcb(Aref(m), m.V)
   dE = zeros(size(Y))
   for el in elements(m.geom.tri)
      if vols[el] > 0
         dE[:, el.t] += vols[el] * stress2frc(el, Wcb1(∇u(el, Y), m.V, Rcb))
      end
   end
   return dE
end

cb_energy1(m::ACModel, Y::Matrix) = cb_energy1(m, Y, m.volT)
