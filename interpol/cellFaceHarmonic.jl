function cellFaceHarmonic(ϕ::CellVariable)
    #=========================================
    Determines the value of the flow variable
    ϕ at the cell faces via harmonic average
    for a uniform mesh
    =========================================#

    dim = ϕ.domain.dimension

    if dim == 1
        dx   = ϕ.domain.cellSize.x
        xval = ϕ.val[1:end-1].*ϕ.val[2:end].*(dx[1:end-1]+dx[2:end]) ./
                (dx[2:end].*ϕ.val[1:end-1]+dx[1:end-1].*ϕ.val[2:end])
        yval = [ ]
        zval = [ ]
    #if dim == 2
    end

    ϕFaceHarmonic = FaceVariable(ϕ.domain, xval, yval, zval)

    return ϕFaceHarmonic
end
