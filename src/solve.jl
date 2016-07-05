
module Solve

using Optim

import LujiaLt: Model, free_defm_indices, defm2dofs, dofs2defm,
            grad, evaluate
import LujiaLt.Potentials: nndist, rdim
import LujiaLt.MDTools: NeighbourList, sites

export solve, quick_solve

"""
`preconditioner(m::Model, Y::Matrix)`
`preconditioner(m::Model, dof::Vector)`

Preconditioner for a standard model (where the dofs are the
   displacements or positions)

* m : model description (e.g., Atm)
* Y : rdim x N matrix of (e.g.) positions
"""
function preconditioner(m::Model, Y::Matrix)
   # get a nearest-neighbour list (+ a bit)
   rnn = nndist(m)
   nlist = NeighbourList(Y, 1.5 * rnn)
   # allocate a triplet format
   I = Int[]; J = Int[]; Z = Float64[];
   # loop through sites and neighbours: add interacting pairs to the matrix
   for (n, neigs, r, R) in sites(nlist), (j, r) in zip(neigs, r)
      z = exp(-3.0*(r/rnn-1.0))
      append!(I, [n;n;j;j]);
      append!(J, [j;n;n;j]);
      append!(Z, [-z;z;-z;z])
   end
   # create sparse matrix and stabilise a little bit
   P = sparse(I, J, Z, size(Y,2), size(Y, 2))
   P += 0.01 * speye(size(Y,2))
   return kron(P, eye(rdim(m)))
end

function preconditioner(m::Model, dof::Vector)
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
               tol = 1e-5,
               display_result = false,
               Optimiser = ConjugateGradient,
               show_trace=false)
   # this needs more boiler plate later on, but for now we can just
   # minimise
   #
   # TODO: depending on the model call a minimise function
   #       or an fzero function
   #
   # TODO: replace m.Yref with reference_configuration ???
   #
   x0 = defm2dofs(m, copy(m.Yref))
   if randomise > 0
      x0 += randomise * 2.0 * (rand(size(x0)) - 0.5)
   end
   obj = x -> evaluate(m, x)
   obj1 = (x, out) -> copy!(out, grad(m, x))
   P = preconditioner(m, x0)
   result = Optim.optimize( DifferentiableFunction(obj, obj1), x0,
                              method = Optimiser(P=P), g_tol=tol,
                              f_tol=0.0, iterations=1_000,
                              show_trace=show_trace )
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




end




"""
`quick_solve(; kwargs...)`

From supplied information, construct model, optimise, and return
equilibrated positions.

## Keyword arguments

* `X` : d x Nat positions matrix
* `Ifree` : list of free atom indices (rest are clamped)
* `Rfree` : ignored if Ifree is provided; otherwise Rfree is the radius
of free atoms
* `V = LennardJonesPotential` : interaction potential
"""
function quick_solve( ; X=nothing, Rfree=nothing, Ifree=nothing,
                        V=Potentials.LennardJonesPotential(),
                        show_trace=false, display_result=true,
                        tol=1e-5 )  # Optimiser = ConjugateGradient,
   # construct free indices (if not supplied)
   if Ifree == nothing
      if Rfree == nothing
         error("`quick_solve` needs either Ifree or Rfree")
      end
      Ifree = find( dist(X) .< Rfree )
   end
   # construct the atomistic model
   model = Atm(X, Ifree, V)
   # solve it
   return Solve.solve(model, show_trace=show_trace,
                     display_result=display_result,
                     tol=tol)  # Optimiser = Optimiser,
end
