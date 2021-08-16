using LinearAlgebra
using Plots

include("../../structs/structs.jl")
include("../../bcs/bcs.jl")
include("../../discret/discret.jl")
include("../../interpol/interpol.jl")
include("../../mesh/mesh.jl")
include("../../solvers/linearSolver.jl")
include("../../navier_stokes/calcCoeff.jl")
include("../../navier_stokes/momentum.jl")
include("../../navier_stokes/RhieChow.jl")

function unsteadyLidDrivenCavity()
    #====================================================
    DESCRIPTION:
    The incompressible Navier-Stokes equations at the
    viscous limit are solved in a 2D boxed domain using
    the SIMPLE algorithm and Rhie-Chow interpolation
    method for a collocated grid.

    ∇⋅(μ∇ū) = ∇p
    ∇⋅ū = 0

    Numerical results are then compared with reference
    data given by Ghia et. al.:
     "High-Re solutions for incompressible flow using the Navier-Stokes
     equations and a multigrid method."
     Journal of computational physics 48.3 (1982): 387-411.

    =====================================================#
##
    #User Input Data
    pRelax      = 1.0                       #Pressure relaxation
    velRelax    = 0.9                       #Velocity relaxation
    uInit       = 0.0                       #Initial guess for x velocity
    vInit       = 0.0                       #Initial guess for y velocity
    pInit       = 0.0                       #Initial guess for pressure
    lidVel      = 1.0                       #Velocity of the lid
    Re          = 100                       #Reynolds number
    cavityLen   = 0.1                       #Length of square cavity
    rho         = 1000                      #Density of fluid
    mu          = rho*cavityLen*lidVel/Re   #Dynamic viscosity
    nIter       = 100                         #Number of iterations
    nx          = 20                        #Number of x nodes
    ny          = 20                        #Number of y nodes
    t           = 100
    dt          = .05

    #Formulate mesh
    msh = meshGen2D(nx, ny, cavityLen, cavityLen)

    #Create and set corresponding velocity boundary conditions
    uBC = generateBC(msh);  vBC = generateBC(msh)

    uBC.top.neu[:] .= 0;    uBC.top.dir[:] .= 1;    uBC.top.val[:] .= lidVel  #Lid velocity
    uBC.bottom.neu[:] .= 0; uBC.bottom.dir[:] .= 1; uBC.bottom.val[:] .= 0    #No slip
    uBC.left.neu[:] .= 0;   uBC.left.dir[:] .= 1;   uBC.left.val[:] .= 0      #No penetration
    uBC.right.neu[:] .= 0;  uBC.right.dir[:] .= 1;  uBC.right.val[:] .= 0     #No penetration

    vBC.top.neu[:] .= 0;    vBC.top.dir[:] .= 1;    vBC.top.val[:] .= 0       #No penetration
    vBC.bottom.neu[:] .= 0; vBC.bottom.dir[:] .= 1; vBC.bottom.val[:] .= 0    #No penetration
    vBC.left.neu[:] .= 0;   vBC.left.dir[:] .= 1;   vBC.left.val[:] .= 0      #No slip
    vBC.right.neu[:] .= 0;  vBC.right.dir[:] .= 1;  vBC.right.val[:] .= 0     #No slip

    #Set the pressure correction equation boundary conditions
    pBC = generateBC(msh)           #Set to Neuman BC by default

    pBC.bottom.neu[div(nx,2)] = 0  #Set reference pressure value
    pBC.bottom.dir[div(nx,2)] = 1
    pBC.bottom.val[div(nx,2)] = 0

    #Build cell and face variables based off initial velocity and pressure guess
    muFaceVar   = generateFaceVar(msh, mu)
    rhoFaceVar  = generateFaceVar(msh, rho)
    faceVel     = generateFaceVar(msh, [uInit vInit])

    uCellVar    = generateCellVar(msh, uInit, uBC)
    vCellVar    = generateCellVar(msh, vInit, vBC)
    pCellVar    = generateCellVar(msh, pInit, pBC)

    #Take the gradient of the pressure
    pGradx, pGrady = grad(pCellVar)

    #main loop
    for j = 1:length(0:dt:t)

        M_ut, RHS_ut = transient(uCellVar, dt)
        M_vt, RHS_vt = transient(vCellVar, dt)
        M_t          = M_ut + M_vt
        RHS_t        = RHS_ut + RHS_vt

        for i = 1:nIter
            #Define previous velocity iteration
            uOld = uCellVar
            vOld = vCellVar

            #Solve momentum equations
            uCellVar, vCellVar, u_ap, v_ap = momentum(msh, uOld, vOld, uBC, vBC,
                                                        pGradx, pGrady, faceVel,
                                                        rho, muFaceVar, velRelax,
                                                        "SIMPLE")

            #Pressure Correction Equation
            # ∇(Vₚ/aₚ)⋅∇p' = ∇⋅ū

            #Rhie Chow Interpolation to find Face Velocity
            faceVel = RhieChow(uCellVar, vCellVar, u_ap, v_ap, pCellVar)

            #Find the RHS value of the pressure correction by taking
            #divergence of the faceVel
            divFaceVel = divergence(faceVel)

            #face values of the velocity coefficients
            u_apFace = arithmeticFaceAvg(u_ap)
            v_apFace = arithmeticFaceAvg(v_ap)

            #Calculate the diffusion coefficient
            diffCoeffx = velRelax*msh.cellSize.y[1]./u_apFace.x
            diffCoeffy = velRelax*msh.cellSize.x[1]./v_apFace.y

            diffCoeff  = FaceVariable(msh, diffCoeffx, diffCoeffy, [])

            #Calculate diffusion term of pressure correction equation
            presCorrectionDiff = diffusionCD(diffCoeff)
            #diff = divergence(grad(diffCoeff))
            #println(presCorrectionDiff.x == diff)

            #Apply pressure boundary conditions
            M_pBC, RHS_pBC = applyBC(pBC)
            # There should be no need to impose pressure BCs to maintain continuity (at least in the fractional step method with a single pressure solve).
            # Instead the corrected face velocities should also satisfy the BCs.
            # Subtracting these from the original face velocities should give the correct boundary velocity corrections and these can be imposed.
            # The pressure field will then be floating with an additive constant, but otherwise the pressure gradients will be correct.

            #Solve for corrected pressure
            pNew = linearSolver(msh, presCorrectionDiff + M_pBC,
                                    divFaceVel[1] + RHS_pBC)

            #Update the old and new pressures and calculate max error
            pOld     = pCellVar
            pCellVar = CellVariable(pCellVar.domain, pCellVar.val + pNew.val*pRelax)

            #Update the velocity given new pressure
            pGradx, pGrady = grad(pNew)

            uNew = uCellVar.val - pGradx.domain.cellSize.y.*pGradx.val./u_ap.val

            vNew = vCellVar.val - pGrady.domain.cellSize.x.*pGrady.val./v_ap.val

            uCellVar = CellVariable(uCellVar.domain,
                        velRelax*uNew + (1-velRelax)*uCellVar.val)
            vCellVar = CellVariable(vCellVar.domain,
                        velRelax*vNew + (1-velRelax)*vCellVar.val)
        end

        #Contour plot of lid driven cavity domain every 50 iterations
        if j%100 == 0
            #pyplot()
            display(contourf(msh.cellCenters.x, msh.cellCenters.y,
                       hypot.(uCellVar.val[2:end-1,2:end-1]',
                       vCellVar.val[2:end-1,2:end-1]'),
                       ylims   = (0,0.1),
                       xlims   = (0,0.1),
                       clims   = (0,1),
                       xlabel  = "x (m)",
                       ylabel  = "y (m)"))

        end


    end
    return uCellVar, vCellVar
end

u, v = unsteadyLidDrivenCavity()
