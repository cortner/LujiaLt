


### This is an attempt to use call-overloading and macros to
### allow the following syntax:
###   pp  ... any subtype of SitePotential or Model
###   p1, p2, ..., pn ... some parameters
###   pp(p1, p2, ..., pn) >>> evaluate(pp, p1, ..., pn)
###   @GRAD pp(r) >>> evaluate(pp, p1, ..., pn, Val{:D}) >>> grad(pp, p1, ..., pn)
###
# create an alias to allow the potential to be called directly
#  make it an inline to ensure that nothing is lost here!
#
# BETTER: call(::Type{Val{:D}}, pp::SimpleFunction, varargs...) = evaluate_d(pp, varargs...)
#
import Base.call
typealias FunUn Union{SitePotential, Model}
@inline call(pp::SitePotential, varargs...) = evaluate(pp, varargs...)
@inline call(pp::SitePotential, ::Type{Val{:GRAD}}, varargs...) = grad(pp, varargs...)

# """`@GRAD` : If `p` is a `SitePotential` or `Model`
# then  `@GRAD p(varargs...)` redirects to `grad(p, varargs...). E.g., a
# basic steepest descent iteration can be written as
# ```
# for n = 1:maxit
#     x = x - step * (@GRAD mymodel(x))
# end
# ```
# """
# macro GRAD(fsig::Expr)
#     try
#         @assert fsig.head == :call
#     catch
#         @show fsig.head
#         @show fsig
#     end
#     for n = 1:length(fsig.args)
#         fsig.args[n] = esc(fsig.args[n])
#     end
#     insert!(fsig.args, 2, Val{:GRAD})
#     return fsig
# end
# export @GRAD
