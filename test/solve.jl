


function test_solve()
   println()
   println("==================================================================")
   println("      TEST_SOLVE")
   println("==================================================================")
   println("Testing implementation of the solver interface")
   model = Atm(V = ToyEAMPotential(), Ra=5.1, defect=:none)
   Ysol = LL.Solve.solve(model, randomise = 0.02, display_result=true)
   return vecnorm(Ysol - model.Yref, Inf) < 1e-4
end




function test_quick_solve()
   println()
   println("==================================================================")
   println("      TEST QUICK_SOLVE")
   println("==================================================================")
   X, _ = LL.lattice_ball(R = 5.1)
   X = X[:, find( dist(X) .> 0.1 )]
   Y = LL.quick_solve(X=X, Rfree=3.1)
   return true
end
