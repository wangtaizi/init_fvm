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
