
module Plotting

using Compose

import LujiaLt.Domain

"return a reasonable choice for the plotting axis"
function autoaxis(X)
    xLim = [extrema(X[1,:])...]
    width = xLim[2]-xLim[1]
    yLim = [extrema(X[2,:])...]
    height = yLim[2]-yLim[1]
    buffer = min( 1.0, 0.05 * width, 0.05 * height )
    xLim[1] -= buffer; xLim[2] += buffer
    yLim[1] -= buffer; yLim[2] += buffer
    return [xLim[1]; xLim[2]; yLim[1]; yLim[2]]
end

_ub_(axis) = UnitBox(axis[1], axis[3], axis[2]-axis[1], axis[4]-axis[3])

"create a context displaying circles for the positions given in `Y`"
function compose_atoms(X; radii=[0.2],
                    fillcolor="tomato", linecolor="darkblue", lwidth=0.3,
                    axis = autoaxis(X) )
    # make sure the radii are an array
    if length(radii) == 1
        radii = [radii...]
    end
    
    # draw the nodes
    ctx = context(units = _ub_(axis) )
    ctx = compose( ctx, circle(X[1,:], X[2,:], radii), 
                   fill(fillcolor), stroke(linecolor), linewidth(lwidth) )
end

    
"create a context displaying the elements (triangles) defined"
function compose_elements(X, T;
                          fillcolor="aliceblue", linecolor="darkblue",
                          lwidth=0.7, axis=autoaxis(X) )
    
    points = Tuple{Float64, Float64}[]
    for n = 1:size(T, 2)
        p = [X[:, T[:, n]] X[:, T[1,n]]]
        for m = 1:size(p, 2)
            push!(points, tuple(p[1, m], p[2,m]))
        end
        push!(points, tuple(NaN, NaN))
    end
    ctx = compose(context(units=_ub_(axis)),
                  polygon(points), 
                  fill(fillcolor), stroke(linecolor), linewidth(lwidth) )
    return ctx
end

    

end
