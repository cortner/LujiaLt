
function test_tb()
   println()
   println("==================================================================")
   println("      TEST TB ASSEMBLY ROUTINES")
   println("==================================================================")
   X, _ = LL.lattice_ball(R = 3.1)
   Y = X + 0.1 * rand(size(X))
   tbm = TBMorsePotential()

   @assert LL.Testing.fdtest( z->tb_energy(tbm, z),
                               z->tb_energy1(tbm, z),
                               Y )
   # site energies test
   println("Testing site energies . . .")
   @assert abs( sum(tb_energies(tbm, Y)) - tb_energy(tbm, Y) ) < 1e-12
   # energy difference test:
   println("Testing energy-difference . . .")
   @assert abs( (tb_energy(tbm, Y) - tb_energy(tbm, X)) -
               tb_energy_diff(tbm, Y, X) ) < 1e-12

   return true
end



function test_tbsolve()
   println()
   println("==================================================================")
   println("      TEST TB SOLVE")
   println("==================================================================")
   println("equilibrate a homogeneous lattice")
   model = TBModel(V=TBMorsePotential(), Rq=5.1, Rbuf=5.1)
   Ysol = LL.Solve.solve(model, randomise = 0.02, display_result=true)
   @assert vecnorm(Ysol - model.Yref, Inf) < 1e-4
   # do a defect equilibration
   println("Now a defect equilibration")
   model = TBModel(V=TBMorsePotential(), Rq=5.1, Rbuf=5.1, defect=:vacancy)
   Ysol = LL.Solve.solve(model, randomise = 0.02, display_result=true)
   return true
end
