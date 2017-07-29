

type BTB
   alpha::Float64
   r0::Float64
   Rcut::Float64
end

BTB(; alpha=2.0, r0=1.1, Rcut=1.8) = BTB(alpha, r0, Rcut)

BTBModel(; alpha=2.0, r0=1.1, Rcut=1.8, beta = 0.1, mu = 0.0 ) =
   TBModel( BTB(alpha, r0, Rcut), mu, beta )


hop(H::BTB, r) = btb_h(r, H.alpha, H.r0, H.Rcut)
hop1(H::BTB, r) = btb_h1(r, H.alpha, H.r0, H.Rcut)
cutoff( H::BTB ) = H.Rcut

## Cutoff function in the NRL ansatz for Hamiltonian matrix elements
## function cutoff(r, Rc, lc)
# Input:
#    r  : variable
#    Rc : cutoff radius
#    lc : cutoff weight
#-------------------------------
# scalar version
btb_cutoff(r::Float64, Rc, lc) = r >= Rc-0.01 ? 0.0 : 1.0 / ( 1.0 + exp( lc / (Rc-r) ) )
# vectorised version
# cutoff(r::Array{Float64}, Rc, lc) = Float64[ cutoff(r[i], Rc, lc) for i in eachindex(r) ]
## first order derivative
#------------------------
btb_cutoff1(r, Rc, lc) = r >= Rc-0.01 ? 0.0 : -1.0 / (1.0 + exp(lc/(Rc-r)))^2 * exp(lc/(Rc-r)) * (lc / (Rc-r)^2)
# vectorised version
# cutoff1(r::Array{Float64}, Rc, lc) = Float64[ cutoff1(r[i], Rc, lc) for i in eachindex(r) ]

## Morse potential
btb_morse(r::Float64, alpha, r0) = exp(-2*alpha*(r-r0)) - 2 * exp(-alpha*(r-r0))
# morse(r::Array{Float64}, alpha, r0) = Float64[ morse(r[i], alpha, r0) for i in eachindex(r) ]
btb_morse1(r::Float64, alpha, r0) = -2*alpha* (exp(-2*alpha*(r-r0)) - exp(-alpha*(r-r0)))
# morse1(r::Array{Float64}, alpha, r0) = Float64[ morse1(r[i], alpha, r0) for i in eachindex(r) ]

## hamiltonian entries potential
btb_h(r, alpha, r0, Rcut) = btb_morse(r, alpha, r0) .* btb_cutoff(r, Rcut, 1.0)
btb_h1(r, alpha, r0, Rcut) = btb_morse1(r, alpha, r0) .* btb_cutoff(r, Rcut, 1.0) + btb_morse(r, alpha, r0) .* btb_cutoff1(r, Rcut, 1.0)
