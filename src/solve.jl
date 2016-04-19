
module Solve


# steepest descent with quadratic line-search
function minimise(;obj=nothing, grad=nothing, x0=nothing,
                  alpha0 = 1.0, tol = 1e-4, precon = 1.0,
                  maxnit = 1_000, Carmijo = 0.2,
                  displevel = 2, alpha_min = 1e-8)

    # evaluate the objective
    x = x0
    E = obj(x)
    ∇E = grad(x)
    Ediff = -0.0
    alpha = 0.0

    if displevel >= 2
        @printf("  nit |    ΔE       |∇E|∞      α \n")
        @printf("------|--------------------------\n")
    end
    for nit = 1:maxnit
        if displevel >= 2
            @printf("%4d | %4.2e  %4.2e  %4.2e \n",
                    nit, Ediff, norm(∇E, Inf), alpha)
        end
        # if residual is less than tolerance >> terminate
        if norm(∇E, Inf) < tol
            if displevel > 0
                println("`minimise` success: residual tolerance reached. ")
            end
            return x, E, ∇E
        end
        
        if nit > 1
            # TODO TODO 
            # extrapolate the step-length (TODO)
            # probably do a secant or BB update?
            alpha = alpha * 1.1
        elseif alpha0 == :auto
            alpha = min(1.0, 1.0/norm(∇E, Inf))
        else
            alpha = alpha0
        end

        # compute energy-difference
        # note that, if `obj` returns site energies, then the
        # sum_kbn([Enwq; -E]) is a robust way to determine the
        # energy-difference. If just total energies are computed,
        # then this is equivalent to `Ediff = Enew - E`.
        ∇E_sq = sumabs2(∇E)
        Enew = obj(x - alpha * ∇E)
        Ediff = sum_kbn([Enew; -E])   
        while Ediff > - alpha * Carmijo * ∇E_sq
            # TODO TODO 
            alpha = alpha / 4
            if alpha < alpha_min
                warn("`minimise` failure: alpha < alpha_min")
                return x, E, ∇E
            end
            Enew = obj(x - alpha * ∇E)
            Ediff = sum_kbn([Enew; -E])   
        end
        
        # do the official update
        x = x - alpha * ∇E
        Eold = E
        E = Enew
        ∇Eold = ∇E
        ∇E = grad(x)
    end

    # the only way to exit from the for-loop is if nit == maxnit
    warn("`minimise` failure: nit > max_nit")
    return x, E, ∇E
end



# CG or BFGS with Armijo + restart



# basic static iteration with force-based (non-)linesearch?



end
