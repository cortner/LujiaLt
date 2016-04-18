
module Potentials

export SitePotential
export evaluate, grad
export ToyEAMPotential
export cutoff
export rDim


"""
Derived types must implement

* evaluate
* grad
"""
abstract SitePotential


########################## CUTOFF FUNCTION #################################

"""A basic C^{2,1} cut-off potential; Arguments:

* r : bond-length (Float or array of floats)
* r0 : inner cut-off radius
* r1 : outer cut-off radius.

`cut` returns `p( (r-r0) / (r1-r0) )`, where

* p(s) = 1, s ≦ 0
* p(s) = 0, s ≧ 1
* p(s) = 6 * (1-s)^5 - 15 * (1-s)^4 + 10 * (1-s)^3

"""
function fcut(r, cutoff)
    s = 1 - (r-cutoff[1]) / (cutoff[2]-cutoff[1])
    return ((s .>= 1) + (0 .<= s .< 1) .* (6 * s.^5 - 15 * s.^4 + 10 * s.^3) )
end

"""Derivative of `cut`; see documentation of `cut`."""
function fcut1(r, cutoff)
    s = 1-(r-cutoff[1]) / (cutoff[2]-cutoff[1])
    return ( - (30*s.^4 - 60 * s.^3 + 30 * s.^2) / (cutoff[2]-cutoff[1])
             .* (0 .< s .< 1) )
end


########################## TOY EAM MODEL #################################

"""
E_n = ∑_ρ ϕ(D_ρ) + F( ∑_ρ ψ(D_ρ) )
ϕ(r) = exp(-2 A (r-1)) - 2 exp(- A(r-1))
ψ(r) = exp(- B (r-1))
F(t) = 0.5 (t-psi0)^2 + C 0.25 (t-psi0)^4
"""
type ToyEAMPotential <: SitePotential
    A
    B
    C
    psi0
    cutoff
end


function ToyEAMPotential(; A=4.0, B=3.0, C=5.0, psi0=6*exp(- B*0.1), cutoff=(1.5, 2.1))
    ToyEAMPotential(A, B, C, psi0, cutoff)
end

cutoff(pot::ToyEAMPotential) = pot.cutoff[2]
rDim(pot::ToyEAMPotential) = 2

morse(r, A) = exp(-2*A*(r-1)) - 2 * exp(-A*(r-1))
morse1(r, A) = -2*A*(exp(-2*A*(r-1)) - exp(-A*(r-1)))
morse(r, A, cutoff) = morse(r, A) .* fcut(r, cutoff)
morse1(r, A, cutoff) =
    morse1(r, A) .* fcut(r, cutoff) + morse(r, A) .* fcut1(r, cutoff)
team_embed(t, C) = 0.5 * t^2 + C * 0.25 * t^4
team_embed1(t, C) = t + C * t^3
team_eldens(r, B, cutoff) = sum(exp(-B*(r-1)) * fcut(r, cutoff))
team_eldens1(r, R, B, cutoff) =
    R .* (exp(-B*(r-1)) .* (fcut1(r, cutoff) - B * fcut(r, cutoff)) ./ r)'

evaluate(p::ToyEAMPotential, r, R) = sum(morse(r, p.A, p.cutoff))

    # (sum(morse(r, p.A, p.cutoff)) +
    #  team_embed( team_eldens(r, p.B, p.cutoff), p.C )
    #  )
    
grad(p::ToyEAMPotential, r, R) =
    R .* (morse1(r, p.A, p.cutoff)./r)'

    # ( R .* (morse1(r, p.A, p.cutoff)./r)' +
    #   team_embed1( team_eldens(r, p.B, p.cutoff), p.C ) *
    #   team_eldens1(r, R, p.B, p.cutoff)
    #   )



end
