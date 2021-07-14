function momentum(msh, uOld, vOld, uBC, vBC, p, faceVel, rho, mu, velRelax, ALGORITHM)
    #===========================================
    DESCRIPTION:
    Given boundary conditions, initial velocity
    components, interpolated face velocity,
    pressure, and algorithm type, solve the
    incompressible navier-stokes momentum
    equation for the updated velocity

    RETURNS:
    CellVariables of updated velocity components
    in addition to their respective coefficients
    ===========================================#

    #Boundary conditions
    M_uBC, RHS_uBC  = applyBC(uBC)
    M_vBC, RHS_vBC  = applyBC(vBC)

    #Discretize diffusive  and convective terms`
    M_diff  = diffusionCD(mu)
    M_conv  = rho*upwindConvection(faceVel)

    #Build total LHS matrix
    u_M = -M_diff + M_conv + M_uBC
    v_M = -M_diff + M_conv + M_vBC

    #Extract coefficients from matrix
    u_ap = calcCoef(msh, u_M, ALGORITHM)
    v_ap = calcCoef(msh, v_M, ALGORITHM)

    #Create cell variables for the pressure
    #gradient in each direction
    pGradx, pGrady = grad(p)

    #Build the RHS vectors for u and v
    u_RHS = RHS_uBC - sourceTerm(pGradx) + (1-velRelax)./velRelax.*diag(u_M).*
                reshape(uOld.val, :, 1)

    v_RHS = RHS_vBC - sourceTerm(pGrady) + (1-velRelax)./velRelax.*diag(v_M).*
                reshape(vOld.val, :, 1)

    #Return the diagonal elements with the velocity
    #relaxation factor implemented
    u_M[diagind(u_M)] = diag(u_M)/velRelax
    v_M[diagind(v_M)] = diag(v_M)/velRelax

    #Solve for the velocity components
    u = linearSolver(msh, u_M, u_RHS)
    v = linearSolver(msh, v_M, v_RHS)

    return u, v, u_ap, v_ap
end
