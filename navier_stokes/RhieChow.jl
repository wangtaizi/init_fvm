function RhieChow(u::CellVariable, v::CellVariable, u_ap::CellVariable,
                    v_ap::CellVariable, p::CellVariable)
    #Returns face values of given field variable
    #within a 2D cartestian grid using the collated
    #Rhie-Chow interpolation method

    #Interpolated face average of the velocity components
    uFace       = arithmeticFaceAvg(u)
    #print(uFace.x[1:5])
    vFace       = arithmeticFaceAvg(v)

    #Interpolated face average of the ap coefficients of
    #the velocity values
    u_apFace     = arithmeticFaceAvg(u_ap)
    v_apFace     = arithmeticFaceAvg(v_ap)

    #cell and face gradients of pressure
    presFaceGrad                = faceGrad(p)
    presFaceGrad.x[:,[1 end]]  .= 0
    presFaceGrad.y[[1 end],:]  .= 0
    #println(presFaceGrad.x)

    presGradx, presGrady        = grad(p)
    #println(presGradx.val)

    #interpolated face average of the cell pressure gradient
    presGradxFace                = arithmeticFaceAvg(presGradx)
    presGradxFace.x[:,[1 end]]  .= 0

    presGradyFace                = arithmeticFaceAvg(presGrady)
    presGradyFace.y[[1 end], :] .= 0

    #Calculate cell volumes and take arithmetic average
    vol     = generateCellVar(p.domain,
                    p.domain.cellSize.x * p.domain.cellSize.y')
    volFace = arithmeticFaceAvg(vol)

    #Define and update starred face velocities
    faceVelx        = uFace.x + volFace.x ./ u_apFace.x .*
                    (presGradxFace.x - presFaceGrad.x)
    faceVelx[1,:]  .= 0
    faceVelx[end,:].= 0

    faceVely        = vFace.y + volFace.y ./ v_apFace.y .*
                    (presGradyFace.y - presFaceGrad.y)
    faceVely[:,1]  .= 0
    faceVely[:,end].= 0

    faceVel    = FaceVariable(u.domain, faceVelx, faceVely, [])


    return faceVel
end
