function generateNodeVar(meshVar::MeshStructure, nodeVal)
    #=================================================
    DESCRIPTION:
    Generate a node variable as a function of
    given mesh and domain properties

    RETURNS:
    NodeVariable with uniform value across domain
    ==================================================#
    rank    = meshVar.dimension
    dim     = meshVar.dims


    if rank == 1
        val  = nodeVal.*ones(dim[1]+1,1)

    else rank == 2
        val  = nodeVal.*ones(dim[1]+1,dim[2]+1)
    end

    nodeVar = NodeVariable(meshVar, val)

    return nodeVar
end
