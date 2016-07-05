

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
   Ysol = LL.Solve.solve(model, display_result=true)
   #  LL.Solve.steepest_descent( obj = z->evaluate(model, z),
   #                     grad = z->grad(model, z),
   #                     x0 = x, alpha0 = 1e-2, maxnit = 100)
   println("==================================================================")
   return passed
end




function test_bqce()
   println()
   println("==================================================================")
   println("      TEST_BQCE")
   println("==================================================================")
   println("Testing implementation of the BQCE model")
   model = BQCE(V = ToyEAMPotential(), Ra=1.1, Rb=2.1, Rc=6.1)
   println(" nX = ", length(find(model.volX .> 0)),
            "; nT = ", length(find(model.volT .> 0)) )
   Y = model.Yref + 0.01 * rand(size(model.Yref))
   println("E = ", evaluate(model, Y))
   println("|dE|_∞ = ", norm( grad(model, Y), Inf ) )
   println("check visually that the finite-difference test looks reasonable.")
   passed = LL.Testing.fdtest( z->evaluate(model, z),
                               z->grad(model, z),
                               Y )
   # println("trying to optimise the geometry")
   # Ysol = LL.Solve.solve(model, display_result=true)
   println("==================================================================")
   return true
end
