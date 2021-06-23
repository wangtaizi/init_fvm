function ghostCells(ϕ::CellVariable, BC::BoundaryCondition)
    #==========================================
    Given the cell values of some parameter ϕ,
    determine the ghost cells based on boundary
    conditions of domain
    ==========================================#

    dim = BC.domain.dimension

    if dim == 1
        ghostϕ = ghostCells_1D(ϕ, BC)

    #if dim == 2
    #    ghostϕ = ghostCells_2D(ϕ, BC)
    end

    return ghostϕ
end

function ghostCells_1D(ϕ::CellVariable, BC::BoundaryCondition)
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
