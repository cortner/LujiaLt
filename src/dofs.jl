
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
