
import .Potentials: SitePotential, cutoff, evaluate, grad
import .MDTools: NeighbourList, sites

# export Model, Atm


########################## atomistic assembly #################################

"""
* `V::SitePotential` : site energy encoding the atomistic model
* `X::Matrix{Float64}` : positions
* `vol` : effective volumes (default 1)
"""
function at_energy(V::SitePotential, X::Matrix{Float64}, vol::Vector{Float64})
    nlist = NeighbourList(X, cutoff(V))
    E = zeros(size(X, 2))
    for (n, neigs, r, R) in sites(nlist)
        if vol[n] <= 0.0; continue; end    # skip vol-0 sites
        E[n] = vol[n] * V(r, R)
    end
    return sum_kbn(E)
end

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
        dV = grad(V, r, R)
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
