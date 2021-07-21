function linearSolver(meshStruct::MeshStructure, M, RHS)
    #====================================================
    Solve a linear system using Julia's built in \ method
    ====================================================#

    n       = meshStruct.dimension
    dims    = meshStruct.dims

    sol     = M \ RHS

        if n >= 2
            ϕ_val   = reshape(sol, dims.+2)
        else
            ϕ_val   = reshape(sol, (dims[1]+2, 1))
        end

    ϕ_struct = CellVariable(meshStruct, ϕ_val)

    return ϕ_struct
end
