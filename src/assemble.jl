
import .Potentials: SitePotential, cutoff, evaluate, grad
import .MDTools: NeighbourList, sites

# export Model, Atm


########################## atomistic assembly #################################

"""
* `V::SitePotential` : site energy encoding the atomistic model
* `X::Matrix{Float64}` : positions
* `vol` : effective volumes (default 1)
"""
function at_energy(V::SitePotential, X::Matrix{Float64}; vol=ones(size(X,2)))
    nlist = NeighbourList(X, cutoff(V))
    E = 0.0
    for (n, neigs, r, R) in sites(nlist)
        if vol[n] <= 0.0; continue; end    # skip vol-0 sites
        E += vol[n] * evaluate(V, r, R)
    end
    return E
end

"""
* `V::SitePotential` : site energy encoding the atomistic model
* `X::Matrix{Float64}` : positions
* `vol` : effective volumes (default 1)
"""
function at_energy1(V::SitePotential, X)
    nlist = NeighbourList(X, cutoff(V))
    dE = zeros(X)
    for (n, neigs, r, R) in sites(nlist)
        dV = grad(V, r, R)
        dE[:, neigs] += dV
        dE[:, n] -= sum(dV, 2)
    end
    return dE
end



######################

# abstract Model

# type Atm <: Model
#     geom::Domain
#     V::SitePotential
#     Ifree::Vector{Int}
#     vol::Vector{Float64}
# end

# function Atm(geom::Domain, V::SitePotential)
    
# end


# dof_vector(m::Atm) = zeros(rDim(m.V) * length(m.Ifree))

# energy(m::Atm, dofs::Vector{Float64}) = at_energy(m.V, dofstodefm(dofs); vol=m.vol)
