


module TightBinding

using ..Potentials
import Potentials: fcut, fcut1, swcut, swcut1

export TBHamiltonian

abstract TBHamiltonian

type

    
## hamiltonian entries potential
h(r, alpha, r0, Rcut) = morse(r, alpha, r0) .* cutoff(r, Rcut, 1.0)
h1(r, alpha, r0, Rcut) = morse1(r, alpha, r0) .* cutoff(r, Rcut, 1.0) + morse(r, alpha, r0) .* cutoff1(r, Rcut, 1.0)

# hamiltonian
# ============
# Compute the objects needed for the BabyTB model with Hamiltonian matrix
#    H_ij = exp(-alpha*rij), i \neq j
#         = 0, i = j
# and repulsive pair potential
#    phi(r) = e0 * exp(- gamma r )
# Parameters:
#       x : positions, d x N array
#   tasks : may be { "H", "dH", "hH", "P", "dP" }
#   alpha, gamma, e0 : model parameters; cf above.
#
function hamiltonian( x, tasks; alpha = 2.0, r0 = 1.0, Rcut = 1.8 )
    
    # extract dimensions
    d, N = size(x)
    # compute generic arrays distances
    R, S = TB.distance_matrices(x, 0)
    
    # assemble whatever is requested
    # first write tasks into a list (if it isn't already)
    if (typeof(tasks) == ASCIIString) || (typeof(tasks)==Char)
        tasks = (tasks,)
    end
    
    # and create a returns list
    ret = ()
    for id = 1:length(tasks)
        
        ## H : HAMILTONIAN
        if tasks[id] == "H"
            H = Float64[ h(S[n,m], alpha, r0, Rcut) * sign(S[n,m]) for n = 1:N, m = 1:N ]
            ret = tuple(ret..., H)
            
        ## dH : Derivative of hamiltonian
        elseif tasks[id] == "dH"
            dH = Float64[ h1(S[n,m], alpha, r0, Rcut) * R[a,n,m] * (sign(S[n,m]) / (S[n,m]+eps()))
                          for a = 1:d, n=1:N, m=1:N ]
            ret = tuple(ret...,dH)

        ## hH : Second derivative of hamiltonian
        elseif tasks[id] == "hH"    #  NOT TESTED
            error("BTB: hH not implemented")
            
        ## P : Pair potential
        elseif tasks[id] == "P"
            ret = tuple(ret..., zeros(N, N))
            
        ## dP : Derivative of pair potential
        elseif tasks[id] == "dP"
            # NOTE: remember that dP_:ij = force of Pij bond acting onto i
#             dEp = - gamma * e0 * Ep
#             dP = [ dEp[m,n]*R[a,m,n]/(S[m,n]+eps())
#                   for a = 1:d, m = 1:N, n = 1:N ]
            dP = zeros(d, N, N)
            ret = tuple(ret...,dP)
            
        else 
            throw(ArgumentError("BTB.hamiltonian: illegal `tasks` argument"))
        end

    end # for id = 1:length(tasks)

    # return the constructed tuple
    # (or its first and only element, if only one is requested)
    if length(tasks) == 1
        return ret[1]
    else
        return ret
    end
    
end # function hamiltonian


end
