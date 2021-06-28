function sourceTerm(ϕ::CellVariable)
    dim = ϕ.domain.dimension

    if  dim == 1
        source = sourceTerm_1D(ϕ::CellVariable)
    elseif dim == 2
        source = sourceTerm_2D(ϕ::CellVariable)

    return source
end

function sourceTerm_1D(ϕ::CellVariable)
    #Returns vector for given constant source term
    #for a 2D domain
end

function sourceTerm_2D(ϕ::CellVariable)
    #Returns vector for given constant source term
    #for a 2D domain

    #Import data from mesh
    nx      = ϕ.domain.dims[1]
    ny      = ϕ.domain.dims[2]
    nodes   = reshape(1:(nx+2)*(ny+2), (nx+2,ny+2))

    #initialize source term vector
    source = zeros((nx+2)*(ny+2), 1)

    i       = reshape(nodes[2:nx+1,2:ny+1], (nx*ny,1))
    source[i] = reshape(ϕ.val[2:nx+1,2:ny+1], (nx*ny,1))

    return source
end
