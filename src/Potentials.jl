
module Potentials

export SitePotential, StandardSitePotential, PairPotential
export ToyEAMPotential, LennardJonesPotential
# export fcut, fcut1, swcut, swcut1
export cutoff

import LujiaLt: evaluate, grad, rdim, nndist


"""
Derived types must implement

* evaluate
* grad
"""
abstract SitePotential

abstract PairPotential <: SitePotential

# revisit whether this is at all needed!?!?
"site potentials described by standard deformation"
abstract StandardSitePotential <: SitePotential

"site potentials for anti-plane displacements"
abstract AntiplaneSitePotential <: SitePotential

# some auxiliary function to be called if only R is available but not r
evaluate(V::SitePotential, R::Matrix) = evaluate(V, sqrt(sumabs2(R,1)[:]), R)
grad(V::SitePotential, R::Matrix) = grad(V, sqrt(sumabs2(R,1)[:]), R)

########################## CUTOFF FUNCTION #################################

"""
A basic C^{2,1} cut-off potential; Arguments:

* r : bond-length (Float or array of floats)
* r0 : inner cut-off radius
* r1 : outer cut-off radius.

`cut` returns `p( (r-r0) / (r1-r0) )`, where

* p(s) = 1, s ≦ 0
* p(s) = 0, s ≧ 1
* p(s) = 6 * (1-s)^5 - 15 * (1-s)^4 + 10 * (1-s)^3

"""
function fcut(r::Real, cutoff)
    s = 1 - (r-cutoff[1]) / (cutoff[2]-cutoff[1])
    s³  = s * s * s
    f = s³ * (10.0 + s * (-15.0 + s * 6.0))
    return (s >= 1) + (0 <= s < 1) * f
end

fcut(r::AbstractVector, cutoff) = [fcut(t, cutoff) for t in r]


"""
Derivative of `cut`; see documentation of `cut`.
"""
function fcut1(r::Real, cutoff)
    s = 1-(r-cutoff[1]) / (cutoff[2]-cutoff[1])
    f = s*s * (30.0 + s * (-60.0 + s * 30.0))
    return (- f) / (cutoff[2]-cutoff[1]) * (0 .< s .< 1)
end

fcut1(r::AbstractVector, cutoff) = [fcut1(t, cutoff) for t in r]


function fcut2(r, cutoff)
   s = 1-(r-cutoff[1]) / (cutoff[2]-cutoff[1])
   return ((4*30)*s.^3 - (3*60) * s.^2 + (2*30) * s.^1) / (cutoff[2]-cutoff[1])^2 .* (0 .< s .< 1)
end


### Unused so far, but keep it here for future use
# """SW-type cutoff function  `function cutoff(r, Rc, lc)`
#
# ## Input:
# *    r  : variable
# *    Rc : cutoff radius
# *    lc : cutoff weight
# """
# swcut(r::Float64, Rc, lc) =
#     (r >= Rc-0.01) ? 0.0 : 1.0 / ( 1.0 + exp( lc / (Rc-r) ) )
# swcut(r::Array{Float64}, Rc, lc) =
#     Float64[ swcut(r[i], Rc, lc) for i in eachindex(r) ]
# "first derivative of `swcut`"
# swcut1(r, Rc, lc) =
#     ( (r >= Rc-0.01) ? 0.0 :
#       -1.0 / (1.0 + exp(lc/(Rc-r)))^2 * exp(lc/(Rc-r)) * (lc / (Rc-r)^2) )
# swcut1(r::Array{Float64}, Rc, lc) =
#     Float64[ swcut1(r[i], Rc, lc) for i in eachindex(r) ]


########################## Some Defaults (can be over-loaded if need be)
cutoff(pot::SitePotential) = pot.cutoff[2]
rdim(pot::SitePotential) = 2
nndist(pot::SitePotential) = 1.0


########################## Lennard-Jones MODEL ###############################

type LennardJonesPotential <: PairPotential
    cutoff::NTuple{2, Float64}
    σ::Float64
end

LennardJonesPotential(;cutoff=(1.5, 2.1), σ=1.0) = LennardJonesPotential(cutoff, σ)

function lj(r::Real)
   t = 1/r
   t³ = t*t*t
   t⁶ = t³ * t³
   return (t⁶ - 2.0) * t⁶
end

function lj1(r::Real)
   t = 1/r
   t³ = t*t*t
   t⁶ = t³ * t³
   return -12.0 * (t⁶ - 1.0) * t⁶ * t
end

lj(r::AbstractVector) = [lj(s) for s in r]
lj1(r::AbstractVector) = [lj1(s) for s in r]
lj2(r) = (12*13) * r.^(-14) - (7*12) * r.^(-8)
lj(r, cutoff) = lj(r) .* fcut(r, cutoff)
lj1(r, cutoff) = lj1(r) .* fcut(r, cutoff) + lj(r) .* fcut1(r, cutoff)
lj2(r, cutoff) = lj2(r)  .* fcut(r, cutoff) + 2 * lj1(r) .* fcut1(r, cutoff) + lj(r) .* fcut2(r, cutoff)

evaluate(p::LennardJonesPotential, r, R) = sum( lj(r/p.σ, p.cutoff) )
grad(p::LennardJonesPotential, r, R) = (R/p.σ) .* (lj1(r/p.σ, p.cutoff) ./ (r/p.σ))' / p.σ


########################## Morse Potential MODEL ###############################

type MorsePotential <: PairPotential
   A::Float64
   cutoff::NTuple{2, Float64}
end

MorsePotential(;A=4.0, cutoff=(1.5, 2.1)) = MorsePotential(A, cutoff)

morse(r, A) = exp(-2.0*A*(r-1)) - 2.0 * exp(-A*(r-1.0))
morse1(r, A) = -2.0*A*(exp(-2.0*A*(r-1.0)) - exp(-A*(r-1.0)))
morse(r, A, cutoff) = morse(r, A) .* fcut(r, cutoff)
morse1(r, A, cutoff) =
    morse1(r, A) .* fcut(r, cutoff) + morse(r, A) .* fcut1(r, cutoff)

evaluate(p::MorsePotential, r, R) = sum(morse(r, p.A, p.cutoff))
grad(p::MorsePotential, r, R) = R .* (morse1(r, p.A, p.cutoff)./r)'


########################### TOY EAM MODEL ####################################

"""E_n = ∑_ρ ϕ(D_ρ) + F( ∑_ρ ψ(D_ρ) )
ϕ(r) = exp(-2 A (r-1)) - 2 exp(- A(r-1))
ψ(r) = exp(- B (r-1))
F(t) = 0.5 (t-psi0)² + C 0.25 (t-psi0)⁴
"""
type ToyEAMPotential <: StandardSitePotential
    A::Float64
    B::Float64
    C::Float64
    psi0::Float64
    cutoff::NTuple{2, Float64}
end

# default constructor
ToyEAMPotential(; A=4.0, B=3.0, C=1.0, psi0=6*exp(B*0.03), cutoff=(1.5, 2.1)) =
    ToyEAMPotential(A, B, C, psi0, cutoff)

team_embed(t, C, t0) = 0.5 * (t/t0-1.0)^2 + C * 0.25 * (t/t0-1.0)^4
team_embed1(t, C, t0) = ((t/t0-1.0) + C * (t/t0-1.0)^3)/t0
team_eldens(r, B, cutoff) = sum(exp(-B*(r-1)) .* fcut(r, cutoff))
team_eldens1(r, R, B, cutoff) =
    R .* (exp(-B*(r-1)) .* (fcut1(r, cutoff) - B * fcut(r, cutoff)) ./ r)'

evaluate(p::ToyEAMPotential, r, R) = (sum(morse(r, p.A, p.cutoff))  +
               team_embed( team_eldens(r, p.B, p.cutoff), p.C, p.psi0 ) )

grad(p::ToyEAMPotential, r, R) = ( R .* (morse1(r, p.A, p.cutoff)./r)' +
                     team_embed1( team_eldens(r, p.B, p.cutoff), p.C, p.psi0 ) *
                     team_eldens1(r, R, p.B, p.cutoff) )


########################## Gupta Potential #################################
# TODO ???implement the Gupta potential for copper???


end
