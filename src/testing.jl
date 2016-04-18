
module Testing

import ..Potentials.SitePotential
import ..at_energy
import ..at_energy1

"generic finite-difference test for scalar F"
function fdtest(F::Function, dF::Function, x)
    E = F(x)
    dE = dF(x)
    # loop through finite-difference step-lengths
    @printf("---------|----------- \n")
    @printf("    h    | error \n")
    @printf("---------|----------- \n")
    for p = 3:10
        h = 0.1^p
        dEh = zeros(dE)
        for n = 1:length(dE)
            x[n] += h
            dEh[n] = (F(x) - E) / h
            x[n] -= h
        end
        @printf(" %1.1e | %4.2e  \n", h, vecnorm(dE - dEh, Inf))
    end
    @printf("---------|----------- \n")
end

"finite-difference test of a SitePotential V"
fdtest(V::SitePotential, X) =
    fd_test(x->at_energy(V, x), x->at_energy1(V, x), X)


end
