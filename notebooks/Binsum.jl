module Binsum

export binsum

function binsum{TI <: Integer}(p::Vector, i::Vector{TI}, maxi::Integer)
    s = zeros(maxi)
    for n in 1:length(p)
        s[i[n]] += p[n]
    end
    return s
end

binsum{TI <: Integer}(p::Vector, i::Vector{TI}) = binsum(p, i, maximum(i))

end