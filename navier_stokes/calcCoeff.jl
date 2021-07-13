using LinearAlgebra

function calcCoeff(msh, M, ALGORITHM)
    #===========================================
    DESCRIPTION:
    Calculates the sum of coeffient terms aₚ for
    the discretized momentum equation given
    the LHS matrix and algorithm type

    RETURNS:
    CellVariable with matrix of summed coeffient
    terms aₚ
    ===========================================#

    nx  = msh.dims[1]
    ny  = msh.dims[2]

    if lowercase(ALGORITHM) == "simple"
        apSum   = diag(M)
    end

    #Coefficients used are normalized by volume,
    #thus to match literature the coeffients
    #must be multiplied by it
    vol = repeat(msh.cellSize.x, 1, ny+2).*repeat(msh.cellSize.y', nx+2, 1)
    val = reshape(apSum, nx+2, ny+2).*vol

    apSum = CellVariable(msh.domain, val)

    return apSum
end
