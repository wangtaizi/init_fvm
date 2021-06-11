function generateBC(meshStruct::MeshStructure)
    #==================================
    General function for formulation of
    Neuman boundary conditions by
    creating BC structure from
    mesh structure
    ==================================#

    dim = meshStruct.dimension

    if dim == 1
        BC = generateBC_1D(meshStruct)
    #elseif dim == 2
    #    BC = generateBC_2D(meshStruct)
    end

    return BC
end
