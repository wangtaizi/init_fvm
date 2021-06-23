function diffusionCD(D::FaceVariable)
    #=======================================
    Given the diffusion coefficient, returns
    a discretized diffusion term using
    the method of central differences
    =======================================#

    dim = D.domain.dimension

    if d == 1
        M = diffusionCD_1D(D)
    end

    return M
end

function diffusionCD_1D(D::FaceVariable)
    #========================================
    Returns a discretized diffusion term in 1D
    using the method of central differences
    ========================================#

    #Mesh data
    nx  = D.domain.dims[1]
    dom = [1:Nx+2]
    D_x = D.domain.cellSize.x # x data of mesh
    dx  = 0.5*(D_x[1:end-1] + D_x[2:end])
    Dx  = D.xval #x component of face value, way too many dx's here :(

    #initialize vectors used to define sparse matrix
    ii  = zeros(3*(nx+2), 1)
    jj  = zeros(3*(nx+2), 1)
    s   = zeros(3*(nx+2), 1)

    #Solve for left and right face variable components
    #that is advected through each cell
    D_l = Dx[1:nx]./(dx[1:nx].*D_x[2:nx+1])
    D_r = Dx[2:nx+1]./(dx[2:nx+1].*D_x[2:nx+1])

    #Formulate internal cells coefficients
    val_l  = reshape(D_l, [nx,1])
    val_r  = reshape(D_r, [nx,1])
    val_p = -(val_l + val_r)

    #Sparse matrix generation
    ii[1:3*nx] = repeat(reshape(dom[2:nx+1], [nx, 1]), outer = [3, 1])
    jj[1:3*nx] = [reshape(dom[1:nx], [nx, 1]); reshape(dom[2:nx+1], [nx, 1])
                    reshape(dom[3:nx+2], [nx, 1])]
    s[1:3*nx]  = [val_l; val_p; val_r]

    M          = sparse(ii[1:3*nx], jj[1:3*nx], s[1:3*nx], [nx+2, nx+2])

    return M
end
