#==================================================
Define function methods to generate a cell variable
as a function of dimensionality of input.
Initial function method does not include ghost cells,
latter method includes this as long as the BC
is also given
==================================================#


function generateCellVar(meshVar::MeshStructure, cellVal)
    #=====================================================
    Generate cell variable based on dimension of problem
    without inclusion of of ghost cell
    =====================================================#
    BC  = generateBC(meshVar)

    dim = meshVar.dims
    if prod(size(cellVal)) == 1 #Number of cellVal elems = 1
        cb  = ghostCells(cellVal*ones(dim), BC)

    elseif size(cellVal) == dim
        cb  = ghostCells(cellVal, BC)

    elseif size(cellVal) == dim.+2
        cb  = cellVal

    else
        cb  = ghostCells(zeros(dim), BC)
    end

    cellVar = CellVariable(meshVar, cb)

    return cellVar
end

function generateCellVar(meshVar::MeshStructure, cellVal, BC::BoundaryCondition)
    #=====================================================
    Generate cell variable based on dimension of problem
    including ghost cells as a function of input BCs
    =====================================================#
    dim = meshVar.dims

    if prod(size(cellVal)) == 1 #Number of cellVal elems = 1
        cb  = ghostCells(cellVal*ones(dim), BC)

    elseif prod(size(cellVal) == dim)
        cb  = ghostCells(cellVal, BC)

    elseif prod(size(cellVal) == dim+2)
        cb  = cellVal

    else
        cb  = ghostCells(zeros(dim), BC)
    end

    cellVar = CellVariable(meshVar, cb)

    return cellVar
end
