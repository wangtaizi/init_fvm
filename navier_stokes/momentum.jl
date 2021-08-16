function momentum(msh, uOld, vOld, uBC, vBC, pGradx, pGrady, faceVel,
                    rho, muFaceVar, velRelax, ALGORITHM)
    #===========================================
    DESCRIPTION:
    Given boundary conditions, initial velocity
    components, interpolated face velocity,
    pressure, and algorithm type, solve the
    incompressible steady navier-stokes momentum
    equation for the updated velocity

    RETURNS:
    CellVariables of updated velocity components
    in addition to their respective coefficients
    ===========================================#

    #Boundary conditions
    M_uBC, RHS_uBC  = applyBC(uBC)
    M_vBC, RHS_vBC  = applyBC(vBC)

    #Discretize diffusive  and convective terms`
    M_diff  = diffusionCD(muFaceVar)
    M_conv  = rho .* upwindConvection(faceVel)[1] # For the two-phase solver later, rho and mu will also be cell variables

    #Build total LHS matrix
    u_M = -M_diff + M_conv + M_uBC
    v_M = -M_diff + M_conv + M_vBC

    #Extract coefficients from matrix
    u_ap = calcCoeff(msh, u_M, ALGORITHM)
    v_ap = calcCoeff(msh, v_M, ALGORITHM)

    #Create cell variables for the pressure
    #gradient in each direction
    #pGradx, pGrady = grad(pCellVar)

    #Build the RHS vectors for u and v
    u_RHS = RHS_uBC - sourceTerm(pGradx) + (1-velRelax)./velRelax.*diag(u_M).*
                reshape(uOld.val, :, 1)

    # Not sure where the velocity relaxation is coming from. Are you using this to mimic time advancement?
    # I see now, this comes from the SIMPLE algorithm. This is OK when the problem is not stiff (so large time steps are OK) or when the problem is steady.
    # But if many small time steps are required then it is better to reduce the number of pressure Poisson solves per time step and a fractional step method might be better.
    # In the fractional step method the pressure term does not appear in the momentum equation and is completely shifted to the projection step with one Poisson solve per time step. But it becomes explicit (as opposed to implicit in SIMPLE).

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

function momentum(msh, uOld, vOld, uBC, vBC, pGradx, pGrady, faceVel,
                    rho, muFaceVar, velRelax, ALGORITHM, M_t, RHS_t)
    #===========================================
    DESCRIPTION:
    Given boundary conditions, initial velocity
    components, interpolated face velocity,
    pressure, and algorithm type, solve the
    incompressible unsteady navier-stokes momentum
    equation for the updated velocity

    RETURNS:
    CellVariables of updated velocity components
    in addition to their respective coefficients
    ===========================================#

    #Boundary conditions
    M_uBC, RHS_uBC  = applyBC(uBC)
    M_vBC, RHS_vBC  = applyBC(vBC)

    #Discretize diffusive  and convective terms`
    M_diff  = diffusionCD(muFaceVar)
    M_conv  = rho .* upwindConvection(faceVel)[1] # For the two-phase solver later, rho and mu will also be cell variables

    #Build total LHS matrix
    u_M = -M_diff + M_conv + M_uBC + M_t
    v_M = -M_diff + M_conv + M_vBC + RHS_t

    #Extract coefficients from matrix
    u_ap = calcCoeff(msh, u_M, ALGORITHM)
    v_ap = calcCoeff(msh, v_M, ALGORITHM)

    #Create cell variables for the pressure
    #gradient in each direction
    #pGradx, pGrady = grad(pCellVar)

    #Build the RHS vectors for u and v
    u_RHS = RHS_uBC - sourceTerm(pGradx) + (1-velRelax)./velRelax.*diag(u_M).*
                reshape(uOld.val, :, 1)

    # Not sure where the velocity relaxation is coming from. Are you using this to mimic time advancement?
    # I see now, this comes from the SIMPLE algorithm. This is OK when the problem is not stiff (so large time steps are OK) or when the problem is steady.
    # But if many small time steps are required then it is better to reduce the number of pressure Poisson solves per time step and a fractional step method might be better.
    # In the fractional step method the pressure term does not appear in the momentum equation and is completely shifted to the projection step with one Poisson solve per time step. But it becomes explicit (as opposed to implicit in SIMPLE).

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
