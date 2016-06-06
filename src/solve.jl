
module Solve

using Optim

import LujiaLt: Model, free_defm_indices, defm2dofs, dofs2defm,
            grad, evaluate
import LujiaLt.Potentials: nndist, rdim
import LujiaLt.MDTools: NeighbourList, sites

export solve


"""
`preconditioner(m::Model, Y::Matrix)`
`preconditioner(m::Model, dof::Vector)`

Preconditioner for a standard model (where the dofs are the
   displacements or positions)

* m : model description (e.g., Atm)
* Y : rdim x N matrix of (e.g.) positions
"""
@noinline function preconditioner(m::Model, Y::Matrix)
   # get a nearest-neighbour list (+ a bit)
   rnn = nndist(m.V)
   nlist = NeighbourList(Y, 1.5 * rnn)
   # allocate a triplet format
   I = Int[]; J = Int[]; Z = Float64[];
   # loop through sites and neighbours: add interacting pairs to the matrix
   for (n, neigs, r, R) in sites(nlist), (j, r) in zip(neigs, r)
      z = exp(-3.0*(r-rnn))
      push!(I, n); push!(J, j); push!(Z, -z)
      push!(I, n); push!(J, n); push!(Z, z)
      push!(I, j); push!(J, n); push!(Z, -z)
      push!(I, n); push!(J, n); push!(Z, z)
   end
   # create sparse matrix and stabilise a bit
   P = sparse(I, J, Z, size(Y,2), size(Y, 2))
   P += 0.01 * speye(size(Y,2))
   return kron(P, eye(rdim(m.V)))
end

@noinline function preconditioner(m::Model, dof::Vector)
   J = free_defm_indices(m)
   return preconditioner(m, dofs2defm(m, dof))[J, J]
end


"""
`function solve(m::Model; kwargs...)`

## Keyword arguments

* `randomise = 0.0` : randomise the initial condition (mostly for testing)
* `gtol = 1e-6` : gradient tolerance
* `display_result = false` : display solver statistics
* `method = ConjugateGradient` : change the optimiser
"""
function solve(m::Model;
               randomise = 0.0,
               tol = 1e-6,
               display_result = false,
               Optimiser = ConjugateGradient )
   # this needs more boiler plate later on, but for now we can just
   # minimise
   #
   # TODO: depending on the model call a minimise function
   #       of an fzero function
   #
   x0 = defm2dofs(m, copy(m.Yref))
   if randomise > 0
      x0 += randomise * 2.0 * (rand(size(x0)) - 0.5)
   end
   obj = x -> evaluate(m, x)
   obj1 = (x, out) -> copy!(out, grad(m, x))
   P = preconditioner(m, x0)
   result = Optim.optimize( DifferentiableFunction(obj, obj1), x0,
                              method = Optimiser(P=P),
                              ftol=0.0, grtol=tol, iterations=20_000 )
   # TODO: analyse `result` more carefully
   if display_result
      println(result)
   end
   return dofs2defm(m, result.minimum)
end




# function minimise(; obj=nothing, grad=nothing, x0=nothing,
#                   alpha0 = 1.0, tol = 1e-4, precon = 1.0,
#                   maxnit = 1_000, Carmijo = 0.2,
#                   displevel = 2, alpha_min = 1e-8 )
#
#    # interface Optim.jl
#
# end



# function fzero(;obj=nothing, grad=nothing, x0=nothing,
#                   alpha0 = 1.0, tol = 1e-4, precon = 1.0,
#                   maxnit = 1_000, Carmijo = 0.2,
#                   displevel = 2, alpha_min = 1e-8)
#
#    # this will use
# end




# TODO:
#  * move everything to just using Optim.jl
#  * implement CG solver from shewchuck; Matlab code copied below.


# include("solve_old.jl")

end
