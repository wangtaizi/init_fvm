function faceGrad(ϕ::CellVariable)
    #Return ∇ϕ as a rank-1 FaceVariable structure

    dim = ϕ.domain.dimension

    if dim == 1
        grad = faceGrad_1D(ϕ::CellVariable)

    elseif dim == 2
        grad = faceGrad_2D(ϕ::CellVariable)

    end

    return grad
end

function faceGrad_1D(ϕ::CellVariable)
    #===========================================
    DESCRIPTION:
    Calculate 2D gradient of field variable
    from interpolated face values using centered
    difference

    RETURNS:
    FaceVariable of face gradient
    ===========================================#

    nx  = ϕ.domain.dims[1]
    lx  = ϕ.domain.cellSize.x
    dx  = 0.5*(lx[1:end-1] + lx[2:end])
    xval = (ϕ.val[2:nx+2] - ϕ.val[1:nx+1])./dx

    #dx  = ϕFace.domain.faceCenters.x[2:end] -
            #ϕFace.domain.faceCenters.x[1:end-1]
    #xval = (ϕFace.x[3:nx+1] - ϕFace.x[1:nx-1])./
                #(dx[1:end-1] + dx[2:end])
    #xval = (ϕFace.x[3:nx+1] - ϕFace.x[1:nx-1])./
    #            (lx[3:nx+1] - lx[1:nx-1])
    yval = [ ]
    zval = [ ]

    grad = FaceVariable(ϕFace.domain, xval, yval, zval)

    return grad
end

function faceGrad_2D(ϕ::CellVariable)
    #===========================================
    DESCRIPTION:
    Calculate 2D gradient of field variable
    from interpolated face values using centered
    difference

    RETURNS:
    FaceVariable of face gradient
    ===========================================#

    nx  = ϕ.domain.dims[1]
    ny  = ϕ.domain.dims[2]
    lx  = repeat(ϕ.domain.cellSize.x, 1, ny)
    ly  = repeat(ϕ.domain.cellSize.y', nx, 1)
    dx  = 0.5*(lx[1:end-1,:] + lx[2:end,:])
    dy  = 0.5*(ly[:,1:end-1] + ly[:,2:end])

    xval = (ϕ.val[2:nx+2, 2:ny+1] - ϕ.val[1:nx+1,2:ny+1])./dx
    yval = (ϕ.val[2:nx+1, 2:ny+2] - ϕ.val[2:nx+1,1:ny+1])./dy
    zval = [ ]


    grad = FaceVariable(ϕFace.domain, xval, yval, zval)

    return grad
end
