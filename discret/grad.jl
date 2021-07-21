function grad(ϕ::CellVariable)
    #Return ∇ϕ as a rank-1 CellVariable structure

    dim = ϕ.domain.dimension

    if dim == 1
        grad = gradient_1D(ϕ::CellVariable)

        return grad

    elseif dim == 2
        gradx, grady = gradient_2D(ϕ::CellVariable)

        return gradx, grady
    end
end

function gradient_1D(ϕ::CellVariable)
    #===========================================
    DESCRIPTION:
    Calculate 2D gradient of field variable
    from adjascent node values using centered
    differences.

    RETURNS:
    CellVariable for each x and y gradient
    ===========================================#
    nx      = ϕ.domain.dims[1]
    xcells  = ϕ.domain.cellSize.x
    dx      = 0.5*(xcells[1:end-1] + xcells[2:end])

    xval = (ϕ.val[3:nx+2] - ϕ.val[1:nx])./
                (dx[1:end-1]+dx[2:end])

    yval = [ ]

    zval = [ ]

    grad = CellVariable(ϕ.domain, xval)

    return grad

end

function gradient_2D(ϕ::CellVariable)
    #===========================================
    DESCRIPTION:
    Calculate 2D gradient of field variable
    from interpolated node values using centered
    differences.

    RETURNS:
    CellVariable for each x and y gradient
    ===========================================#
    nx      = ϕ.domain.dims[1]
    ny      = ϕ.domain.dims[2]
<<<<<<< HEAD
    xcells  = repeat(ϕ.domain.cellSize.x, 1, ny)
    ycells  = repeat(ϕ.domain.cellSize.y', nx, 1)
=======
    xcells  = repeat(ϕ.domain.cellSize.x, 1,ny)
    ycells  = repeat(ϕ.domain.cellSize.y', nx,1)
>>>>>>> 27d496454d80003064881718e297164123f62914

    dx      = 0.5*(xcells[1:end-1,:] + xcells[2:end,:])
    dy      = 0.5*(ycells[:,1:end-1] + ycells[:,2:end])

    xval = zeros(nx+2, ny+2)
    xval[2:nx+1,2:ny+1] = (ϕ.val[3:nx+2, 2:ny+1] - ϕ.val[1:nx,2:ny+1])./
                (dx[1:end-1,:]+dx[2:end,:])

    yval = zeros(nx+2, ny+2)
    yval[2:nx+1,2:ny+1] = (ϕ.val[2:nx+1, 3:ny+2] - ϕ.val[2:nx+1,1:ny])./
                (dy[:,1:end-1]+dy[:,2:end])

    zval = [ ]

    gradx = CellVariable(ϕ.domain, xval)
    grady = CellVariable(ϕ.domain, yval)

    return gradx, grady
end
