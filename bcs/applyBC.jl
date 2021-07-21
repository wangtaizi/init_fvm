function applyBC(BC::BoundaryCondition)
    #==============================
    Generate matrix of boundary condition
    coefficients and vector of RHS values
    ==============================#

    dim = BC.domain.dimension

    if dim == 1
        bc_mat, bc_rhs = applyBC_1D(BC)
    else dim == 2
        bc_mat, bc_rhs = applyBC_2D(BC)
    end

    return bc_mat, bc_rhs
end

function applyBC_1D(BC::BoundaryCondition)
    #===========================================
    Generate boundary condition matrix of
    coefficients and RHS vector for a 1D domain
    ===========================================#

    #Mesh structure data
    nx      = BC.domain.dims[1]
    dx_i    = BC.domain.cellSize.x[1]
    dx_f    = BC.domain.cellSize.x[end]
    domain  = 1:nx+2
    boundaryNodes = 8

    #Initialize vectors of indicies
    ii  = zeros(boundaryNodes, 1)
    jj  = zeros(boundaryNodes, 1)
    s   = zeros(boundaryNodes, 1)

    #Initalize RHS vector
    bc_rhs = zeros(nx+2, 1)

    #Generate BC matrix and RHS vector

    #Right Boundary Condition
    i       = nx+2
    q       = 0

    q       += 1
    ii[q]   = domain[i]
    jj[q]   = domain[i]
    s[q]    = BC.right.dir/2 + BC.right.neu/dx_f #by the way, the Dirichlet condition here isn't compatible with the harmonic mean in the diffusion test case

    q       += 1
    ii[q]   = domain[i]
    jj[q]   = domain[i-1]
    s[q]    = BC.right.dir/2 - BC.right.neu/dx_f

    bc_rhs[domain[i]] = BC.right.val

    #Left Boundary Condition
    i       = 1

    q       += 1
    ii[q]   = domain[i]
    jj[q]   = domain[i+1]
    s[q]    = -(BC.left.dir/2 + BC.left.neu/dx_i) #left

    q       += 1
    ii[q]   = domain[i]
    jj[q]   = domain[i]
    s[q]    = -(BC.left.dir/2 - BC.left.neu/dx_i) #left

    bc_rhs[domain[i]] = -(BC.left.val)

    #Formulate Sparse Matrix
    bc_matrix = sparse(ii[1:q], jj[1:q], s[1:q], nx+2, nx+2)

    return bc_matrix, bc_rhs
end

function applyBC_2D(BC::BoundaryCondition)
    #===========================================
    Generate boundary condition matrix of
    coefficients and RHS vector for a 1D domain
    ===========================================#

    #Mesh structure data
    nx      = BC.domain.dims[1]
    ny      = BC.domain.dims[2]
    dx_i    = BC.domain.cellSize.x[1]
    dx_f    = BC.domain.cellSize.x[end]
    dy_i    = BC.domain.cellSize.y[1]
    dy_f    = BC.domain.cellSize.y[end]
    domain  = 1:nx+2
    boundaryNodes = 8*(nx+ny+2)
    nodes   = reshape(1:(nx+2)*(ny+2), (nx+2,ny+2))

    #Initialize vectors of indicies
    ii  = zeros(boundaryNodes, 1)
    jj  = zeros(boundaryNodes, 1)
    s   = zeros(boundaryNodes, 1)

    #Initalize RHS vector
    bc_rhs = zeros((nx+2)*(ny+2), 1)

    #Generate BC matrix and RHS vector

    #Assign zeros to unused corner nodes
    q       = 1:4
    ii[q]   = BC.domain.cornerNodes
    jj[q]   = BC.domain.cornerNodes
    s[q]    .= maximum(BC.top.dir/2 + BC.top.neu/dy_f)
    bc_rhs[BC.domain.cornerNodes] .= 0

    #Top Boundary Condition
    i       = 2:nx+1 #i=2:nx+1?
    j       = ny+2 #j=ny+2?
    q       = q[end] .+ (1:nx)
    ii[q]   = nodes[i,j]
    jj[q]   = nodes[i,j]
    s[q]    = BC.top.dir/2 + BC.top.neu/dy_f

    q       = q[end] .+ (1:nx) #q=q[end]+(1:nx)?
    ii[q]   = nodes[i,j]
    jj[q]   = nodes[i,j-1]
    s[q]    = BC.top.dir/2 - BC.top.neu/dy_f

    bc_rhs[nodes[i,j]] = BC.top.val

    #Bottom Boundary Condition
    i       = 2:nx+1 #i=2:nx+1?
    j       = 1 #j=1?
    q       = q[end] .+ (1:nx)
    ii[q]   = nodes[i,j]
    jj[q]   = nodes[i,j+1]
    s[q]    = -(BC.bottom.dir/2 + BC.bottom.neu/dy_i)

    q       = q[end] .+ (1:nx) #q=q[end]+(1:nx)?
    ii[q]   = nodes[i,j]
    jj[q]   = nodes[i,j]
    s[q]    = -(BC.bottom.dir/2 - BC.bottom.neu/dy_i)

    bc_rhs[nodes[i,j]] = -(BC.bottom.val)

    #Right Boundary Condition
    i       = nx+2
    j       = 2:ny+1
    q       = q[end] .+ (1:ny)
    ii[q]   = nodes[i,j]
    jj[q]   = nodes[i,j]
    s[q]    = BC.right.dir/2 + BC.right.neu/dx_f

    q       = q[end] .+ (1:ny)
    ii[q]   = nodes[i,j]
    jj[q]   = nodes[i-1,j]
    s[q]    = BC.right.dir/2 - BC.right.neu/dx_f

    bc_rhs[nodes[i,j]] = BC.right.val

    #Left Boundary Condition
    i       = 1
    j       = 2:ny+1
    q       = q[end] .+ (1:ny)
    ii[q]   = nodes[i,j]
    jj[q]   = nodes[i+1,j]
    s[q]    = -(BC.left.dir/2 + BC.left.neu/dx_i)

    q       = q[end] .+ (1:ny)
    ii[q]   = nodes[i,j]
    jj[q]   = nodes[i,j]
    s[q]    = -(BC.left.dir/2 - BC.left.neu/dx_i)

    bc_rhs[nodes[i,j]] = -(BC.left.val)

    #Formulate Sparse Matrix
    bc_matrix = sparse(ii[1:q[end]], jj[1:q[end]], s[1:q[end]],
                    (nx+2)*(ny+2), (nx+2)*(ny+2))

    return bc_matrix, bc_rhs
end
