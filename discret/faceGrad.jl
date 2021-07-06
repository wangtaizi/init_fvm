function faceGrad(ϕ::CellVariable)
    #Return ∇ϕ as a rank-1 FaceVariable structure

    dim = ϕ.domain.dimension

    if dim == 1
        grad = faceGrad_1D(ϕ::CellVariable)

    elseif dim == 2
        grad = faceGrad_2D(ϕ::CellVariable)

    return grad
end

function faceGrad_1D(ϕ::CellVariable)
    #Calculate the 1D x gradient

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
    ϕFace   = arithmetcFaceAvg(ϕ)

    nx  = ϕFace.domain.dims[1]
    ny  = ϕFace.domain.dims[2]
    lx  = repeat(ϕFace.domain.cellSize.x, 1, ny)
    ly  = repeat(ϕFace.domain.cellSize.y', nx, 1)
    dx  = 0.5*(lx[1:end-1,:] + lx[2:end,:])
    dy  = 0.5*(ly[:,1:end-1] + ly[:,2:end])

    #xval = (ϕ.val[2:nx+2, 2:ny+1] - ϕ.val[1:nx+1,2:ny+1])./dx
    #yval = (ϕ.val[2:nx+1, 2:ny+2] - ϕ.val[2:nx+1,1:ny+1])./dy
    #zval = [ ]

    xval = (ϕFace.x[3:nx+2, 2:ny+1] - ϕFace.x[1:nx,2:ny+1])./
                (dx[1:end-1]+dx[2:end])
    yval = (ϕFace.y[2:nx+1, 3:ny+2] - ϕFace.y[2:nx+1,1:ny])./
                (dy[1:end-1]+dy[2:end])
    zval = [ ]

    grad = FaceVariable(ϕFace.domain, xval, yval, zval)

    return grad
end
