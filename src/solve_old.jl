




# steepest descent with quadratic line-search
function steepest_descent(;obj=nothing, grad=nothing, x0=nothing,
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






fxunction ncg_force(; grad=nothing, x0=nothing, tol = 1e-4,
                     σ_0 = 1.0, maxnf = 1_000, alpha_min = 1e-8,
                     σ0 = 0.1, P = 1.0, precondprep! = (P, x) -> P )

   nf = 0
   k = 0
   x = x0
   r = - grad(x)
   nf += 1
   P = precondprep!(P, x)
   s = similar(r)
   A_ldiv_B!(s, P, r)
   d = copy(s)
   δ_new = dot(r, d)
   δ_0 = δ_new
   res_inf = norm(r, Inf)
   while nf < maxnf && res_inf > tol
      j = 0
      δ_d = dot(d, r)
      α = - σ_0    # σ_0 is undefined >>> input!
      η_prev = dot( grad(x + σ_0 * d), d )



   # allocate some stuff
   N = length(x0)
   u = x0             # current state
   u_aux = copy(x0)

   # Evaluate first force and search direction
   F = grad(u)
   p = - F

   # some initialisation
   δ_new = dot(F, p)
   δ0 = δ_new
   nf = 0
   res_inf = norm(F, Inf)

   while nf < maxnf && res_inf > tol
      # I have no idea what all this is ?!?!?
      sec_iter = 0  # secant iterations?
      δ = dot(p, p)    # delta_D = δ
      α = - σ0
      u_aux = u + σ0 * p
      F_new = grad(u_aux)
      η_prev = dot(F_new, p)
      η = - dot(F, p)
      α *= η / (η_prev - η)
      u += α * p
      η_prev = η
      sec_iter += 1
      #  Perform secant iteration loop
      while sec_iter < max_sec_iter && α^2 * δ > tol^2
#             F = grad(u)
#             η = dot(F, p)
#             α *= η / (η_prev - η)
#             u += α * p
#             η_prev = η
#             sec_iter += 1
#         end



#     end






# #     force = -bqcf(u);

# #     delta_OLD = delta_NEW;
# #     delta_MID = dot(force,sForce);

# #     %In this step would compute a preconditioner
# #     sForce = force;  %would be sForce = M^(-1)force

# #     delta_NEW = dot(force,sForce);
# #     beta = (delta_NEW-delta_MID)/delta_OLD;
# # %     disp('Here!!!!!');
# # %     pause;
# # %     disp('Delta_NEW');
# # %     disp(delta_NEW);
# # %     disp('beta');
# # %     disp(beta);
# # %     pause;
# #     if beta <= 0
# #         direction = sForce;
# #     else
# #         direction = sForce + beta*direction;
# #     end
# #     iter = iter+1;
# #     fprintf('frc : %4.2e,  \n', delta_NEW);
# #     fprintf('delta_D : %4.2e,  \n', delta_D);
# #     fprintf('delta_MID : %4.2e,  \n', delta_MID);
# #     fprintf('delta_OLD : %4.2e,  \n', delta_OLD);
# #     fprintf([' sec_iter = ', num2str(sec_iter)]);
# #     fprintf(' beta = %4.2e,  ', beta);
# #     fprintf(' alpha = %4.2e,  ', alpha);
# #     disp(['|U|_inf = ', num2str(norm(u(:), Inf))]);

# # end

# # disp(iter);
# # disp(delta_NEW);
# # pause;

# # end




# include("solve_old.jl")
