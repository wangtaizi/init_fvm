function arithmeticFaceAvg(ϕ::CellVariable)
    #Return FaceVariable corresponding to arithmetic
    #interpolation of adjacent cells

    dim = ϕ.domain.dimension

    if dim == 1
        dx      = ϕ.domain.cellSize.s
        xval    = (dx[1:end-1].*ϕ.val[1:end-1]+dx[2:end].*ϕ.val[2:end])./
                    (dx[2:end]+dx[1:end-1])
        yval    = [ ]
        zval    = [ ]

    elseif dim == 2
        nx      = ϕ.domain.dims[1]
        ny      = ϕ.domain.dims[2]
        dx      = repeat(ϕ.domain.cellSize.x, 1, ny)
        dy      = repeat(ϕ.domain.cellSize.y', nx, 1)

        xval    = (dy[1:end-1,:].*ϕ.val[1:end-1,2:end-1]+
                    dx[2:end,:].*ϕ.val[2:end,2:end-1])./
                    (dx[2:end,:]+dx[1:end-1,:])

        yval    = (dy[:,1:end-1].*ϕ.val[2:end-1,1:end-1]+
                    dy[:,2:end].*ϕ.val[2:end-1,2:end])./
                    (dy[:,2:end]+dy[:,1:end-1])

        zval    = [ ]
    end

    avg = FaceVariable(ϕ.domain, xval, yval, zval)

    return avg
end
