
"""
`extend(Xlge, X, U)`

extend an approximate solution to be defined
on the nodes of a larger domain defined through m
"""
function extend(Xlge, X, U)
    tri = FEM.Triangulation(X)
    idx = FEM.locate(Xlge, tri)
    Uext = zeros(size(U,1), length(idx))
    for n = 1:length(idx)
        if idx[n] == 0
            Uext[:, n] = 0.0
        else
            t = tri.T[:, idx[n]]
            lam = FEM.convex_coordinates(Xlge[:,n], idx[n], tri)
            Uext[:, n] = U[:, t] * lam
        end
    end
    return Uext
end

fdiff(U, E) = U[:, E[1,:][:]] - U[:, E[2,:][:]]

function H1seminorm(tri, U)
   E = FEM.edges(tri)
   DU = fdiff(U, E)
   DX = fdiff(tri.X, E)
   return (sumabs2(DU, 1) ./ sumabs2(DX, 1)) |> sum |> sqrt
end

"""
`compute_error(m, m_ex, U_ex, E_ex)`

compute the Ḣ¹ error and energy error of the model `m`,
as measured against a computation on a much larger domain
(m_ex, U_ex, E_ex).
"""
function error_energynorm(Y, at, Ylge, atlge)
   X, Xlge = positions(at), positions(atlge)
   U, Ulge = Y - X, Ylge - Xlge
   Uext = extend(Xlge, X, U)
   # compute and return errors
   return H1seminorm(atlge.geom.tri, Ulge - Uext)
end


function error_energynorm(Y, X::Matrix, Ylge, Xlge::Matrix)
   U, Ulge = Y - X, Ylge - Xlge
   Uext = extend(Xlge, X, U)
   tri = FEM.Triangulation(Xlge)
   # compute and return errors
   return H1seminorm(tri, Ulge - Uext)
end
