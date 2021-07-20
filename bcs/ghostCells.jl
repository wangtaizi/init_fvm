function ghostCells(ϕ::Any, BC::BoundaryCondition)
    #==========================================
    Given the cell values of some parameter ϕ,
    determine the ghost cells based on boundary
    conditions of domain
    ==========================================#

    dim = BC.domain.dimension

    if dim == 1
        ghostϕ = ghostCells_1D(ϕ, BC)

    elseif dim == 2
        ghostϕ = ghostCells_2D(ϕ, BC)
    end

    return ghostϕ
end

function ghostCells_1D(ϕ::Any, BC::BoundaryCondition)
    #============================================
    Determine the value of the ghost cells for a 1D
    domain and add them to 1D ϕ domain
    ============================================#

    dx_i = BC.domain.cellSize.x[1]
    dx_f = BC.domain.cellSize.x[end]

    ghostϕ = [(BC.left.val-ϕ[1]*(BC.left.neu/dx_i + BC.left.dir/2)) /
                (-BC.left.neu/dx_i + BC.left.dir/2); ϕ;
                (BC.right.val-ϕ[end]*(-BC.right.neu/dx_f + BC.right.dir/2)) /
                (BC.right.neu/dx_f + BC.right.dir/2)]

    return ghostϕ
end

function ghostCells_2D(ϕ::Any, BC:: BoundaryCondition)
    #============================================
    Determine the value of the ghost cells for a 2D
    domain and add them to 2D ϕ domain
    ============================================#

    #Mesh data
    nx   = BC.domain.dims[1]
    dx_i = BC.domain.cellSize.x[1]
    dx_f = BC.domain.cellSize.x[end]

    ny   = BC.domain.dims[2]
    dy_i = BC.domain.cellSize.y[1]
    dy_f = BC.domain.cellSize.y[end]

    #Initialize output ϕ
    ghostϕ  = zeros(nx+2, ny+2)
    ghostϕ[2:nx+1, 2:ny+1]  = ϕ

    #left boundary
    ghostϕ[1,2:ny+1] = (BC.left.val' .- ϕ[1,:].*(BC.left.neu'/dx_i + BC.left.dir'/2)) ./
                (-BC.left.neu'/dx_i + BC.left.dir'/2)

    #right boundary
    ghostϕ[nx+2,2:ny+1] = (BC.right.val'.-ϕ[end,:].*(-BC.right.neu'/dx_f + BC.right.dir'/2)) ./
                (BC.right.neu'/dx_f + BC.right.dir'/2)

    #top boundary
    ghostϕ[2:nx+1,ny+2] = (BC.top.val.-ϕ[:,end].*(-BC.top.neu/dy_f + BC.top.dir/2)) ./
                (BC.top.neu/dy_f + BC.top.dir/2)

    #bottom boundary
    ghostϕ[2:nx+1,1] = (BC.bottom.val.-ϕ[:,1].*(BC.bottom.neu/dy_i + BC.bottom.dir/2)) ./
                (-BC.bottom.neu/dy_i + BC.bottom.dir/2)

    return ghostϕ
end
