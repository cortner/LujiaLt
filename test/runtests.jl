



using LujiaLt; LL=LujiaLt
using LujiaLt.Potentials




function test_toyeam()
    println()
    println("==================================================================")
    println("      TEST_TOYEAM")
    println("==================================================================")
    println("Testing implementation of ToyEAMPotential, at_energy, at_energy1")
    println("check visually that the finite-difference test looks reasonable.")
    V = ToyEAMPotential(B=3.0, C=1.0)
    geom = Domain(Ra=4.0, defect=:vacancy)
    # perturbed configuration
    X = positions(geom)
    Y = X + 0.1 * rand(X)
    passed = LL.Testing.fdtest(V, Y)
    println("==================================================================")
    return passed
end


function test_Atm()
   println()
   println("==================================================================")
   println("      TEST_ATM")
   println("==================================================================")
   println("Testing implementation of the Atm model")
   model = Atm(V = ToyEAMPotential(), Ra=5.0, defect=:vacancy)
   x = LL.defm2dofs(model, model.Yref);
   println("E = ", evaluate(model, x))
   println("|dE|_∞ = ", norm( grad(model, x), Inf ))
   println("check visually that the finite-difference test looks reasonable.")
   passed = LL.Testing.fdtest( z->evaluate(model, z),
                               z->grad(model, z),
                               x + 0.03 * rand(x) )
   println("trying to optimise the geometry")
   Ysol = LL.Solve.solve(model)
   #  LL.Solve.steepest_descent( obj = z->evaluate(model, z),
   #                     grad = z->grad(model, z),
   #                     x0 = x, alpha0 = 1e-2, maxnit = 100)
   println("==================================================================")
   return passed
end


function test_solve()
   println()
   println("==================================================================")
   println("      TEST_SOLVE")
   println("==================================================================")
   println("Testing implementation of the solver interface")
   model = Atm(V = ToyEAMPotential(), Ra=5.0, defect=:none)
   Ysol = LL.Solve.solve(model, randomise = 0.01, display_result=true)

   # x = LL.defm2dofs(model, model.Yref)
   # x += 0.02 * (rand(size(x)) - 0.5)
   # LL.Solve.steepest_descent( obj = z->evaluate(model, z),
   #                     grad = z->grad(model, z),
   #                     x0 = x, alpha0 = 1e-2, maxnit = 100)

   return vecnorm(Ysol - model.Yref, Inf) < 1e-4
end

#############################################################################

# @assert test_toyeam()
# @assert test_Atm()
@assert test_solve()
