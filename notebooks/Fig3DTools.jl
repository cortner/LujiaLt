
module Fig3DTools


get_r(X::Matrix) = sqrt(sumabs2(X,1))


function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repmat(vx, m, 1), repmat(vy, 1, n))
end

function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T},
                     vz::AbstractVector{T})
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end


## function dgrid(vxyz, d)
# generate d-dimensional grid, return is d-tuple
function dgrid(vxyz, d)
    if d == 1
        return vxyz
    elseif d == 2
        return meshgrid(vxyz, vxyz)
    elseif d == 3
        return meshgrid(vxyz, vxyz, vxyz)
    else
        throw(ArgumentError("dgrid: d must be 1, 2, or 3."))
    end
end

## function dgrid_list
# like dgrid, but returns a (d x npoints) array
function dgrid_list(vxyz, d)
    if d == 2
        order = [2, 1]
    elseif d == 3
        order = [2,1,3]
    end
    grids = dgrid(vxyz, d)
    if d == 1; grids = (grids,); end
    nX = length(grids[1])
    x = zeros(d, nX)
    for α = 1:d
        x[α, :] = reshape(grids[order[α]], (1, nX))
    end
    return x
end



# a standard constructor for the LatticeGeom type
# returns a LatticeGeom with
function lattice_ball(A, R)
    # dimension of the lattice
    dim = size(A, 1)
    # get smallest singular value of A know how big the
    # lattice-box needs to be.
    sig = minimum(svd(A)[2])
    # generate a cubic portion of Zd containing -R/sig : R/sig in each
    # coordinate direction
    ndim = ceil(Int, R/sig)

    # obtain the deformed point list
    Z = dgrid_list(-ndim:ndim, dim)
    # the reference configuration
    X = A * Z
    # find points inside the ball and return them
    I_ball = find( sumabs2(X, 1) .<= R^2 )
    # create the array of nuclei
    X = X[:, I_ball]
    return X
end




A_tri = [ [1.,0.] [cos(pi/3), sin(pi/3)] ]
A_fcc = [ [0.,1.,1.] [1.,0.,1.] [1.,1.,0.] ] / sqrt(2.)
A_bcc = [ [1.,-1.,-1.] [1.,1.,-1.] [1.,1.,1.] ] / sqrt(3.)




"""
Generate an approximate ball in a Bravais lattice, possibly with a
point defect.

INPUT (all are keyword arguments)

* lattice : string, so far have implemented "tri", "fcc", "bcc"
* R : scalar, radius
* defect : string, so far have implemented "", "vac", "int",

OUTPUT: y

* y : dim x nX positions array
"""
function lattice_ball(;lattice="tri", R=10, defect="", A = [])
    # determine orientation matrix
    # and the interstitial position
    if lattice == "tri"
        A = A_tri
        xi = 0.5 * A[:,1]
    elseif lattice == "fcc"
        A = A_fcc
        xi = 0.5 * (A[:,1]+A[:,2]+A[:,3])
    elseif lattice == "bcc"
        A = A_bcc
    elseif A == []
        error("LjGeometry.lattice_ball: unknown lattice")
    end
    # problem dimension
    d = size(A)[1]
    # get smallest singular value
    sig = minimum(svd(A)[2])
    # generate a cubic portion of Zd containing -R/sig : R/sig in each
    # coordinate direction
    ndim = ceil(R/sig)
    z = dgrid_list(-ndim:ndim, d)
    # transform via A
    x = A*z
    # find points inside the ball and return them
    I_ball = find( sumabs2(x, 1) .< R^2 )
    y = x[:, I_ball]
    I0 = find( sumabs(y, 1) .== 0.0 )[1]

    # introduce the defect
    if defect == ""     # no defect
        return y, I0
    elseif defect == "vac"    # vacancy defect
        y =[y[:, 1:(I0-1)] y[:, (I0+1):end]]
        return y
    elseif defect == "int"    # interstitial
        y = [y xi]
        return y
    else
        error("LjGeometry.lattice_ball: unknown defect")
    end
end


function multilattice_ball(;lattice="hex", R=10, A = [], p = [])

    if lattice == "hex"
        l = sin(2*pi/3) / sin(pi/6)
        A = [ [l, 0.] l*[.5,sqrt(3/4)] ]
        p = [0.,1.]
    elseif A == [] || p == []
        error("Unkown lattice")
    end

    y1, I0 = lattice_ball(lattice="x", A = A, R=R+1)
    y2 = y1 .+ p
    y = [y1 y2]
    I_ball = find( sumabs2(y, 1) .< R^2 )

    return y[:, I_ball]
end




"""
Writes the atomic positions stored in Y (dim x n-atoms) into an `filename`, in
the xyz format. The extension `.xyz` is not automatically added.

INPUT

* Y : dim x n Float64 array, each column is the position of an atom
* filename : string
* comment : string for the comment line in the xyz format (optional)

REMARKS

* at the moment there is almost no error handling
"""
function write_xyz(Y::Array{Float64, 2}, filename;
                   comment="no comment", extra=[],
                   species=[])

    # analyze data
    dim, num_atoms = size(Y)
    if (dim < 2) || (dim > 3)
        error("LjIO.write_xyz requires 2D or 3D objects")
    end
    # append zeros if Y is a two-dimensional array
    if dim == 2
        Y = [Y; zeros(1, num_atoms)]
    end

    if isempty(species)
        species = ["No" for n = 1:size(Y, 2)]
    end

    # open file for writing
    file = open(filename, "w")

    # write header
    write(file, string(num_atoms, "\n"))
    write(file, string(comment, "\n"))

    # write the atom positions
    for n = 1:num_atoms
        write(file,
              @sprintf("%s \t %f \t %f \t %f ", species[n], Y[1, n], Y[2,n], Y[3,n])
              )
        if !isempty(extra)
            for i = 1:size(extra, 1)
                write(file, @sprintf("\t %f ", extra[i, n]))
            end
        end
        write(file, "\n")
    end

    # close the file
    close(file)
end


end
