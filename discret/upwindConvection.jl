function upwindConvection(u::FaceVariable)
    #Discretization of the convective term
    #∇⋅(uϕ) using upwind scheme

    dim = u.domain.dimension

    if  dim == 1
        M = upwindConvection_1D(u::FaceVariable)

        return M
    elseif dim == 2
        M, Mx, My = upwindConvection_2D(u::FaceVariable)

        return M, Mx, My

end

function upwindConvection_1D(u::FaceVariable)
    #Discretize convection term using
    #upwind scheme of 1D domain
end

function upwindConvection_2D(u::FaceVariable)
    #Discretize convection term using
    #upwind scheme of 2D cartesian domain

    #Import data from mesh
    nx      = u.domain.dims[1]
    ny      = u.domain.dims[2]
    nodes   = reshape(1:(nx+2)*(ny+2), (nx+2,ny+2))
    u_xp  = repeat(u.domain.cellSize.x[2:end-1], 1, ny)
    u_yp  = repeat(u.domain.cellSize.y[2:end-1]', nx, 1)

    #initialize indexing vectors
    ix = zeros(3*(nx+2)*(ny+2), 1)
    jx = zeros(3*(nx+2)*(ny+2), 1)
    sx = zeros(3*(nx+2)*(ny+2), 1)

    iy = zeros(3*(nx+2)*(ny+2), 1)
    jy = zeros(3*(nx+2)*(ny+2), 1)
    sy = zeros(3*(nx+2)*(ny+2), 1)

    #initialize directional max/min vectors
    ue_min  = u.x[2:nx+1,:]
    ue_max  = u.x[2:nx+1,:]
    uw_min  = u.x[1:nx,:]
    uw_max  = u.x[1:nx,:]
    un_min  = u.y[:,2:ny+1]
    un_max  = u.y[:,2:ny+1]
    us_min  = u.y[:,1:ny]
    us_max  = u.y[:,1:ny]

    ue_min[u.x[2:nx+1,:] .> 0.0]    .= 0.0
    ue_max[u.x[2:nx+1,:] .< 0.0]    .= 0.0
    uw_min[u.x[1:nx,:] .> 0.0]      .= 0.0
    uw_max[u.x[1:nx,:] .< 0.0]      .= 0.0
    un_min[u.y[2:ny+1,:] .> 0.0]    .= 0.0
    un_max[u.y[2:ny+1,:] .< 0.0]    .= 0.0
    us_min[u.y[1:ny,:] .> 0.0]      .= 0.0
    us_max[u.y[1:ny,:] .< 0.0]      .= 0.0

    #Internal cell coefficient values
    ae  = ue_min./u_xp
    aw  = -uw_max./u_xp
    an  = un_min./u_yp
    as  = -us_max./u_yp
    apx = (ue_max-uw_min)./u_xp
    apy = (un_max-us_min)./u_yp

    #Correct for cell values at boundries
    #left boundary
    apx[1,:]    = apx[1,:] - uw_max[1,:]/(2*u_xp[1])
    aw[1,:]     = aw[1,:]/2

    #right
    apx[end,:]  = apx[end,:] + ue_min[end,:]/(2*u_xp[end])
    ae[end,:]   = ae[end,:]/2

    #bottom
    apy[:,1]    = apy[:,1] - us_max[:,1]/(2*u_yp[1])
    as[:,1]     = as[:,1]/2

    #top
    apy[:,end]  = apy[:,end] + un_min[:,end]/(2*u_yp[end])
    an[:,end]   = an[:,end]/2

    #Indexing for developing the sparse matrix
    ix[1:3*(nx*ny)] = repeat(reshape(nodes[2:nx+1,2:ny+1], (nx*ny,1)), 3, 1)
    jx[1:3*(nx*ny)] = [reshape(nodes[1:nx,2:ny+1], (nx*ny,1));
                        reshape(nodes[2:nx+1,2:ny+1], (nx*ny,1));
                        reshape(nodes[3:nx+2,2:ny+1], (nx*ny,1));]
    sx[1:3*(nx*ny)] = [reshape(aw, (nx*ny,1));
                        reshape(apx, (nx*ny,1));
                        reshape(ae, (nx*ny,1))]

    iy[1:3*(nx*ny)] = repeat(reshape(nodes[2:nx+1,2:ny+1], (nx*ny,1)), 3, 1)
    jy[1:3*(nx*ny)] = [reshape(nodes[2:nx+1,1:ny], (nx*ny,1));
                        reshape(nodes[2:nx+1,2:ny+1], (nx*ny,1));
                        reshape(nodes[2:nx+1,3:ny+2], (nx*ny,1));]
    sy[1:3*(nx*ny)] = [reshape(as, (nx*ny,1));
                        reshape(apy, (nx*ny,1));
                        reshape(an, (nx*ny,1))]

    #Define sparse matricies
    Mx  = sparse(ix[1:3(nx*ny)], jx[1:3*(nx*ny)], sx[1:3*(nx*ny)],
            (nx+2)*(ny+2), (nx+2)*(ny+2))

    My  = sparse(iy[1:3(nx*ny)], jy[1:3*(nx*ny)], sy[1:3*(nx*ny)],
            (nx+2)*(ny+2), (nx+2)*(ny+2))

    M   = Mx + My
    
    return M, Mx, My
end
