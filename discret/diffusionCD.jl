using SparseArrays

function diffusionCD(D::FaceVariable)
    #=======================================
    Given the diffusion coefficient, returns
    a discretized diffusion term using
    the method of central differences
    =======================================#

    dim = D.domain.dimension

    if dim == 1
        M = diffusionCD_1D(D)
    elseif dim == 2
        M = diffusionCD_2D(D)
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
    dom = 1:nx+2
    D_x = D.domain.cellSize.x # x data of mesh
    dx  = 0.5*(D_x[1:end-1] + D_x[2:end])
    Dx  = D.x #x component of face value, way too many dx's here :(

    #initialize vectors used to define sparse matrix
    ii  = zeros(3*(nx+2), 1)
    jj  = zeros(3*(nx+2), 1)
    s   = zeros(3*(nx+2), 1)

    #Solve for left and right face variable components
    #that is advected through each cell
    D_l = Dx[1:nx]./(dx[1:nx].*D_x[2:nx+1])
    D_r = Dx[2:nx+1]./(dx[2:nx+1].*D_x[2:nx+1])

    #Formulate internal cells coefficients
    val_l  = reshape(D_l, (nx,1))
    val_r  = reshape(D_r, (nx,1))
    val_p = -(val_l + val_r)

    #Sparse matrix generation
    ii[1:3*nx] = [reshape(dom[2:nx+1], (nx, 1)); reshape(dom[2:nx+1], (nx, 1));
                    reshape(dom[2:nx+1], (nx, 1))]
    jj[1:3*nx] = [reshape(dom[1:nx], (nx, 1)); reshape(dom[2:nx+1], (nx, 1))
                    reshape(dom[3:nx+2], (nx, 1))]
    s[1:3*nx]  = [val_l; val_p; val_r]

    M          = sparse(ii[1:3*nx], jj[1:3*nx], s[1:3*nx], nx+2, nx+2)

    return M
end

function diffusionCD_2D(D::FaceVariable)
    #========================================
    Returns a discretized diffusion term in 2D
    using the method of central differences
    for a uniform grid
    ========================================#

    #Mesh data
    nx  = D.domain.dims[1]
    ny  = D.domain.dims[2]
    nodes   = reshape(1:(nx+2)*(ny+2), (nx+2,ny+2))

    D_x = repeat(D.domain.cellSize.x, 1, ny) # x data of mesh
    D_y = repeat(D.domain.cellSize.y', nx, 1)
    dx  = 0.5*(D_x[1:end-1,:] + D_x[2:end,:])
    dy  = 0.5*(D_y[:,1:end-1] + D_y[:,2:end])
    Dx  = D.x #x component of face value, way too many dx's here :(
    Dy  = D.y

    #initialize vectors used to define sparse matrix
    iix  = zeros(3*(nx+2)*(ny+2), 1)
    jjx  = zeros(3*(nx+2)*(ny+2), 1)
    sx   = zeros(3*(nx+2)*(ny+2), 1)

    iiy  = zeros(3*(nx+2)*(ny+2), 1)
    jjy  = zeros(3*(nx+2)*(ny+2), 1)
    sy   = zeros(3*(nx+2)*(ny+2), 1)

    #Reassign directional velocity vectors
    # of each direction
    D_l = Dx[1:nx,:]./(dx[1:nx,:].*D_x[2:nx+1,:]) #left
    D_r = Dx[2:nx+1,:]./(dx[2:nx+1,:].*D_x[2:nx+1,:]) #right
    D_t = Dy[:,2:ny+1]./(dy[:,2:ny+1].*D_y[:,2:ny+1]) #top
    D_b = Dy[:,1:ny]./(dy[:,1:ny].*D_y[:,2:ny+1]) #bottom

    #Formulate internal cells coefficients
    val_l  = reshape(D_l, (nx*ny,1))
    val_r  = reshape(D_r, (nx*ny,1))
    val_t  = reshape(D_t, (nx*ny,1))
    val_b  = reshape(D_b, (nx*ny,1))

    val_px = reshape(-(val_l + val_r), (nx*ny,1))
    val_py = reshape(-(val_t + val_b), (nx*ny,1))

    #Sparse matrix generation
    iix[1:3*(nx*ny)] = [reshape(nodes[2:nx+1, 2:ny+1], (nx*ny, 1));
                            reshape(nodes[2:nx+1, 2:ny+1], (nx*ny, 1));
                            reshape(nodes[2:nx+1, 2:ny+1], (nx*ny, 1))]
    jjx[1:3*(nx*ny)] = [reshape(nodes[1:nx, 2:ny+1], (nx*ny, 1));
                            reshape(nodes[2:nx+1, 2:ny+1], (nx*ny, 1));
                            reshape(nodes[3:nx+2, 2:ny+1], (nx*ny, 1));]
    sx[1:3*(nx*ny)]  = [val_l; val_px; val_r]


    iiy[1:3*(nx*ny)] = [reshape(nodes[2:nx+1, 2:ny+1], (nx*ny, 1));
                            reshape(nodes[2:nx+1, 2:ny+1], (nx*ny, 1));
                            reshape(nodes[2:nx+1, 2:ny+1], (nx*ny, 1))]
    jjy[1:3*(nx*ny)] = [reshape(nodes[2:nx+1, 1:ny], (nx*ny, 1));
                            reshape(nodes[2:nx+1, 2:ny+1], (nx*ny, 1));
                            reshape(nodes[2:nx+1, 3:ny+2], (nx*ny, 1));]
    sy[1:3*(nx*ny)]  = [val_b; val_py; val_t]

    Mx              = sparse(iix[1:3*(nx*ny)], jjx[1:3*(nx*ny)], sx[1:3*(nx*ny)],
                    (nx+2)*(ny+2), (nx+2)*(ny+2))
    My              = sparse(iiy[1:3*(nx*ny)], jjy[1:3*(nx*ny)], sy[1:3*(nx*ny)],
                                    (nx+2)*(ny+2), (nx+2)*(ny+2))
    M               = Mx + My

    return M
end
