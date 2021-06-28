function gradient(ϕ::CellVariable)
    #Given a field variable, return the gradient

    dim = ϕ.domain.dimension

    if dim == 1
        grad = gradient_1D(ϕ::CellVariable)

    elseif dim == 2
        grad = gradient_2D(ϕ::CellVariable)

    return grad
end

function gradient_1D(ϕ::CellVariable)
    #Calculate the 1D x gradient

end

function gradient_2D(ϕ::CellVariable)
    #calculate 2d gradient in x and y directions
    nx      = ϕ.domain.dims[1]
    ny      = ϕ.domain.dims[2]
    xcells  = repeat(ϕ.domain.cellSize.x, (1,ny))
    ycells  = repeat(ϕ.domain.cellSize.y', (nx,1))
    dx      = 0.5*(xcells[1:end-1,:] + xcells[2:end,:])
    dy      = 0.5*(ycells[:,1:end-1] + ycells[:,2:end])

    xval = (ϕ.val[2:nx+2, 2:ny+1] - ϕ.val[1:nx+1,2:ny+1])./dx
    yval = (ϕ.val[2:nx+1, 2:ny+2] - ϕ.val[2:nx+1,1:ny+1])./dy
    zval = [ ]

    grad = FaceVariable(ϕ.domain, xval, yval, zval)

    return grad
end
