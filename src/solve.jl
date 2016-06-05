
module Solve

using Optim

import LujiaLt: Model, free_defm_indices, defm2dofs, dofs2defm,
            grad, evaluate

"""
`preconditioner(m::Model, Y::Matrix)`
`preconditioner(m::Model, dof::Vector)`

Preconditioner for a standard model (where the dofs are the
   displacements or positions)

* m : model description (e.g., Atm)
* Y : rdim x N matrix of (e.g.) positions
"""
function preconditioner_scalar(m::Model, Y::Matrix)
   # get a nearest-neighbour list (+ a bit)
   rnn = nndist(m.V)
   nlist = NeighbourList(X, 1.5 * rnn)
   # allocate a triplet format
   I = Int[]; J = Int[]; Z = Float64[];
   # loop through sites and neighbours
   for (n, neigs, r, R) in sites(nlist)
      for j in neigs
         push!(I, n); push!(J, j); push!( Z, exp(-3.0*(r-rnn)) )
      end
   end
   # create sparse matrix and stabilise a bit
   P = sparse(I, J, Z, size(Y,2), size(Y, 2)) + 0.01 * speye(size(Y,2))
   return kron(P, eye(rdim(m.V)))
end

function preconditioner(m::Model, dof::Vector)
   P = preconditioner(m, dofs2defm(m, dof))[m.Ifree, m.Ifree]
   J = free_defm_indices(m)
   return P[J, J]
end


"""
`function solve(m::Model; kwargs...)`

## Keyword arguments

* `randomise = 0.0` : randomise the initial condition (mostly for testing)
* `gtol = 1e-6` : gradient tolerance
"""
function solve(m::Model;
               randomise = 0.0,
               tol = 1e-6,
               display_result = false
               )
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
   P = nothing; # preconditioner(m, x0)
   result = Optim.optimize( DifferentiableFunction(obj, obj1), x0,
                              method = ConjugateGradient(P=P),
                              ftol=1e-32, grtol=tol )
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



function fzero(;obj=nothing, grad=nothing, x0=nothing,
                  alpha0 = 1.0, tol = 1e-4, precon = 1.0,
                  maxnit = 1_000, Carmijo = 0.2,
                  displevel = 2, alpha_min = 1e-8)

   # this will use
end




# TODO:
#  * move everything to just using Optim.jl
#  * implement CG solver from shewchuck; Matlab code copied below.


include("solve_old.jl")



end
