function arithmeticNodeAvg(ϕ::CellVariable)
    #Return NodeVariable corresponding to arithmetic
    #interpolation of adjacent cells (linear and bilinear interpolation)

    dim = ϕ.domain.dimension

    if dim == 1
        dx      = ϕ.domain.cellSize.x
        val    = (dx[2:end].*ϕ.val[1:end-1]+dx[1:end-1].*ϕ.val[2:end])./
                    (dx[2:end]+dx[1:end-1])

    elseif dim == 2
        nx      = ϕ.domain.dims[1]
        ny      = ϕ.domain.dims[2]
        dx      = repeat(ϕ.domain.cellSize.x, 1, ny+2)
        dy      = repeat(ϕ.domain.cellSize.y', nx+2, 1)

		val		= (dx[2:end,2:end].*dy[2:end,2:end].*ϕ.val[1:end-1,1:end-1]+
					dx[1:end-1,2:end].*dy[1:end-1,2:end].*ϕ.val[2:end,1:end-1]+
					dx[2:end,1:end-1].*dy[2:end,1:end-1].*ϕ.val[1:end-1,2:end]+
					dx[1:end-1,1:end-1].*dy[1:end-1,1:end-1].*ϕ.val[2:end,2:end])./
					(dx[1:end-1,1:end-1].*dy[1:end-1,1:end-1]+
					dx[2:end,1:end-1].*dy[2:end,1:end-1]+
					dx[1:end-1,2:end].*dy[1:end-1,2:end]+
					dx[2:end,2:end].*dy[2:end,2:end])

    end

    avg = NodeVariable(ϕ.domain, val)

    return avg
end
