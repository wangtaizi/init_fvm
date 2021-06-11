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
