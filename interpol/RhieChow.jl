function RhieChow(u::CellVariable, v::CellVariable, uap::CellVariable,
                    vap::CellVariable, p::CellVariable)
    #Returns face values of given field variable
    #within a 2D cartestian grid using the collated
    #Rhie-Chow interpolation method

    #Interpolated face average of the velocity components
    uFace       = arithmeticAvg(u)
    vFace       = arithmeticAvg(v)

    #Interpolated face average of the ap coefficients of
    #the velocity values
    uapFace     = arithmeticAvg(uap)
    vapFace     = arithmeticAvg(vap)

    #cell and face gradients of pressure
    presFaceGrad            = faceGrad(p)
    presGradx, presGrady    = grad(p)

    #interpolated face average of the cell pressure gradient
    presGradxFace   = arithmeticAvg(presGradx)
    presGradyFace   = arithmeticAvg(presGrady)

    #Define and update starred face velocities
    faceVel    = FaceVariable(u.domain, xval, yval)

    faceVel.x  = uFace.x + (p.domain.cellSize.x .* p.domain.cellSize.y) ./
                    uapFace.x .* (presGradxFace - presFaceGrad.x)

    faceVel.y  = vFace.y + (p.domain.cellSize.x .* p.domain.cellSize.y) ./
                    vapFace.y .* (presGradyFace - presFaceGrad.y)


    return faceVel
end
