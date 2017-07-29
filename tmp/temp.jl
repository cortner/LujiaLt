"""
Compute vector and scalar pairwise distances. This is a faster version
(ca factor 2) of
```
d, N = size(X)
R = Float64[ X[a, i] - X[a,j] for a=1:d, i=1:N, j=1:N ]
S = sqrt(squeeze(sumabs2(R,1), 1))
```
"""
function distance_matrices(X)
   d, N = size(X)
   @assert d == 2
   R = zeros(d, N, N); S = zeros(N, N)
   t = (-1)::Int; s = 0::Int
   @inbounds for j = 1:2:2*N-1, i = 1:2:2*N-1
         s = s+1; t = t+2
         R[t] = X[i]-X[j]
         R[t+1] = X[i+1]-X[j+1]
         S[s] = sqrt(R[t]*R[t] + R[t+1]*R[t+1])
      end
   end
   return R, S
end
 
