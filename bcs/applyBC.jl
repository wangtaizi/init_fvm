function applyBC(BC::BoundaryCondition)
    #==============================
    Generate matrix of boundary condition
    coefficients and vector of RHS values
    ==============================#

    dim = BC.domain.dimension

    if dim == 1
        [bc_mat, bc_rhs] = applyBC_1D(BC)
    #else dim == 2
    #    [bc_mat, bc_rhs] = applyBC_2D(BC)
    end

    return [bc_mat, bc_rhs]
end
