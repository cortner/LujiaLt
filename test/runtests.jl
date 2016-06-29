



using LujiaLt; LL=LujiaLt
using LujiaLt.Potentials
using LujiaLt.TightBinding




function test_toyeam()
    println()
    println("==================================================================")
    println("      TEST_TOYEAM")
    println("==================================================================")
    println("Testing implementation of ToyEAMPotential, at_energy, at_energy1")
    println("check visually that the finite-difference test looks reasonable.")
    V = ToyEAMPotential(B=3.0, C=1.0)
    geom = Domain(Ra=4.1, defect=:vacancy)
    # perturbed configuration
    X = positions(geom)
    Y = X + 0.1 * rand(X)
    passed = LL.Testing.fdtest(V, Y)
    println("==================================================================")
    return passed
end


function test_lennardjones()
    println()
    println("==================================================================")
    println("      TEST_TOYEAM")
    println("==================================================================")
    println("Testing implementation of ToyEAMPotential, at_energy, at_energy1")
    println("check visually that the finite-difference test looks reasonable.")
    V = LennardJonesPotential()
    geom = Domain(Ra=4.1, defect=:vacancy)
    # perturbed configuration
    Y = positions(geom); Y += 0.1 * rand(size(Y))
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
   Ysol = LL.Solve.solve(model, display_result=true)
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
   model = Atm(V = ToyEAMPotential(), Ra=5.1, defect=:none)
   Ysol = LL.Solve.solve(model, randomise = 0.02, display_result=true)
   return vecnorm(Ysol - model.Yref, Inf) < 1e-4
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


import LujiaLt.TightBinding: hop, hop1


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






#############################################################################


# @assert test_toyeam()
# @assert test_lennardjones()
# @assert test_Atm()
# @assert test_solve()
# @assert test_bqce()
# @assert test_quick_solve()
# @assert test_tb()
@assert test_tbsolve()
