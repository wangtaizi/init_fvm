function generateFaceVar(meshVar::MeshStructure, faceVal)
    #=================================================
    DESCRIPTION:
    Generate a face variable as a function of
    given mesh and domain properties

    RETURNS:
    FaceVariable with uniform value across domain
    ==================================================#
    rank    = meshVar.dimension
    dim     = meshVar.dims


    if rank == 1
        xval  = faceVal[1].*ones(dim[1]+1,1)
        yval  = [ ]
        zval  = [ ]

    else rank == 2
        if length(faceVal) == 2
            xval  = faceVal[1].*ones(dim[1]+1,dim[2])
            yval  = faceVal[2].*ones(dim[1],dim[2]+1)
            zval  = [ ]

        else
            xval  = faceVal[1].*ones(dim[1]+1,dim[2])
            yval  = faceVal[1].*ones(dim[1],dim[2]+1)
            zval  = [ ]
        end
    end

    faceVar = FaceVariable(meshVar, xval, yval, zval)

    return faceVar
end
