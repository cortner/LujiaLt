
module Plotting

using Compose

import ..Domain

"return a reasonable choice for the plotting axis"
function autoaxis(X)
    xLim = [extrema(tri.X[1,:])...]
    width = xLim[2]-xLim[1]
    yLim = [extrema(tri.X[2,:])...]
    height = yLim[2]-yLim[1]
    buffer = min( 1.0, 0.05 * width, 0.05 * height )
    xLim[1] -= buffer; xLim[2] += buffer
    yLim[1] -= buffer; yLim[2] += buffer
    return [xLim[1]; xLim[2]; yLim[1]; yLim[2]]
end

"create a context displaying circles for the positions given in `Y`"
function compose_atoms(Y; radii=[0.2],
                       fillcolor="tomato", linecolor="darkblue", lwidth=0.0,
                       axis = autoaxis(Y) )
    # make sure the radii are an array
    if length(radii == 1)
        radii = [radii...]
    end
    
    # draw the nodes
    ctx = context(units = [axis[1]; axis[3]; axis[2]-axis[1]; axis[4]-axis[3]])
    ctx = compose(ctx,
                  circle(tri.X[1,:], tri.X[2,:], radii), 
                  fill(fillcolor), stroke(linecolor) )
end

"create a context displaying the elements (triangles) defined"
function compose_elements()
    
end


end
