function applyBC(BC::BoundaryCondition)
    #==============================
    Generate matrix of boundary condition
    coefficients and vector of RHS values
    ==============================#

    dim = BC.domain.dimension

    if dim == 1
        bc_mat, bc_rhs = applyBC_1D(BC)
    #else dim == 2
    #    [bc_mat, bc_rhs] = applyBC_2D(BC)
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
    domain  = [1:nx+2]
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

    q       += 1
    ii[q]   = domain[i]
    jj[q]   = domain[i]
    s[q]    = BC.right.dir/2 + BC.right.neu/dx_f

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
    s[q]    = -(BC.right.dir/2 + BC.right.neu/dx_i)

    q       += 1
    ii[q]   = domain[i]
    jj[q]   = domain[i-1]
    s[q]    = -(BC.right.dir/2 - BC.right.neu/dx_i)

    bc_rhs[domain[i]] = -(BC.left.val)

    #Formulate Sparse Matrix
    bc_matrix = sparse(ii[1:q], jj[1:q], s[1:q], [nx+2,nx+2])

    return bc_matrix, bc_rhs
end
