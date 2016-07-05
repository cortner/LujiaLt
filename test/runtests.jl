

function test_frchybrid_linear()
   model = FrcHybridModel(Vqm=TBMorsePotential(), Rqm=3.1, Rmm=10.1, Rbuf=5.1)
   tbm = TBModel(model.Xref, model.V, model.Ifree, model.Yref, [0.0])

   println("Testing that the linear MM force is O(h^2) accurate")
   U1 = rand(size(model.Yref))
   err = [];
   for h in [0.1; 0.03; 0.01; 0.003; 0.001]
      Y = model.Yref + h * U1
      fqm = grad(tbm, Y)
      fmm = copy(fqm)
      fmm[:, model.Imm] = grad(model, Y)[:, model.Imm]
      push!(err, vecnorm(fqm-fmm)/h^2)
      @show vecnorm(fqm-fmm)/h^2
   end
   # this asserts that the error doesnt grow as h decreases
   @assert maximum(err) < 2 * err[1]

   return true
end







using LujiaLt; LL=LujiaLt
using LujiaLt.Potentials
using LujiaLt.TightBinding


include("potentials.jl")
include("atc.jl")
include("tb.jl")
include("solve.jl")



#############################################################################


# @assert test_toyeam()
# @assert test_lennardjones()
# @assert test_Atm()
# @assert test_solve()
# @assert test_bqce()
# @assert test_quick_solve()
# @assert test_tb()
# @assert test_tbsolve()
@assert test_frchybrid_linear()


# @assert test_frchybrid_quad()
