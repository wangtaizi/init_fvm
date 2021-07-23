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

function couetteFlow()
    #====================================================
    DESCRIPTION:
    The incompressible Navier-Stokes equations are solved
    for a fully developed cacnonical couette flow
    within a 2D channel
    =====================================================#
##
    #User Input Data
    pRelax      = 1.0                       #Pressure relaxation
    velRelax    = 0.9                       #Velocity relaxation
    uInit       = 0.0                       #Initial guess for x velocity
    vInit       = 0.0                       #Initial guess for y velocity
    pInit       = 0.0                       #Initial guess for pressure
    surfaceVel  = 1.0                       #Velocity of the top surface
    height      = 1.0                       #channel height
    len         = 3*height                  #channel length
    Re          = 100                       #Reynolds number
    rho         = 1000                      #Density of fluid
    mu          = rho*len*surfaceVel/Re     #Dynamic viscosity
    nIter       = 500                       #Number of iterations
    nx          = 100                       #Number of x nodes
    ny          = 100                       #Number of y nodes

    #Formulate mesh
    msh = meshGen2D(nx, ny, height, len)

    #Create and set corresponding velocity boundary conditions
    uBC = generateBC(msh);  vBC = generateBC(msh)

    uBC.top.neu[:] .= 0;    uBC.top.dir[:] .= 1;    uBC.top.val[:] .= surfaceVel  #surface velocity
    uBC.bottom.neu[:] .= 0; uBC.bottom.dir[:] .= 1; uBC.bottom.val[:] .= 0    #No slip
    uBC.left.neu[:] .= 0;   uBC.left.dir[:] .= 1;   uBC.left.val[:] .= 0      #No penetration
    uBC.right.neu[:] .= 0;  uBC.right.dir[:] .= 1;  uBC.right.val[:] .= 0     #No penetration

    vBC.top.neu[:] .= 0;    vBC.top.dir[:] .= 1;    vBC.top.val[:] .= 0       #No penetration
    vBC.bottom.neu[:] .= 0; vBC.bottom.dir[:] .= 1; vBC.bottom.val[:] .= 0    #No penetration
    vBC.left.neu[:] .= 0;   vBC.left.dir[:] .= 1;   vBC.left.val[:] .= 0      #No slip
    vBC.right.neu[:] .= 0;  vBC.right.dir[:] .= 1;  vBC.right.val[:] .= 0     #No slip

    #Set the pressure correction equation boundary conditions
    pBC = generateBC(msh)           #Set to Neuman BC by default

    pBC.top.neu[:] .= 0;    pBC.top.dir[:] .= 1;    pBC.top.val[:] .= 0       #No penetration
    pBC.bottom.neu[:] .= 0; pBC.bottom.dir[:] .= 1; pBC.bottom.val[:] .= 0    #No penetration
    pBC.left.neu[:] .= 0;   pBC.left.dir[:] .= 1;   pBC.left.val[:] .= 0      #No slip
    pBC.right.neu[:] .= 0;  pBC.right.dir[:] .= 1;  pBC.right.val[:] .= 0     #No slip

    #Build cell and face variables based off initial velocity and pressure guess
    muFaceVar   = generateFaceVar(msh, mu)
    rhoFaceVar  = generateFaceVar(msh, rho)
    faceVel     = generateFaceVar(msh, [uInit vInit])

    uCellVar    = generateCellVar(msh, uInit, uBC)
    vCellVar    = generateCellVar(msh, vInit, vBC)
    pCellVar    = generateCellVar(msh, pInit, pBC)

    #main loop
    for i = 1:nIter
        #Define previous velocity iteration
        uOld = uCellVar
        vOld = vCellVar

        #Solve momentum equations
        uCellVar, vCellVar, u_ap, v_ap = momentum(msh, uOld, vOld, uBC, vBC, pCellVar, faceVel,
                                    rho, muFaceVar, velRelax, "SIMPLE")

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
        #println(diag(presCorrectionDiff))

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

        #================================================================#

        #Calculating Error Relative to Ghia Solution
        #Interpolate datapoints
        uitp     = interpolate((msh.cellCenters.y,), uCellVar.val[div(nx,2), 2:end-1],
                        Gridded(Linear()))
        uetp     = extrapolate(uitp, Flat())
        uGhia    = uetp(GhiaValues.y*cavityLen)

        vitp     = interpolate((msh.cellCenters.x,), vCellVar.val[div(ny,2), 2:end-1],
                        Gridded(Linear()))
        vetp     = extrapolate(uitp, Flat())
        vGhia    = vetp(GhiaValues.y*cavityLen)

        if i%50 == 0
            #Calculate root mean squared error
            uRMSE   = sqrt(sum((GhiaValues.u[:,1]-uGhia[:]).^2)/length(GhiaValues.u[:,1]))
            vRMSE   = sqrt(sum((GhiaValues.v[:,1]-vGhia[:]).^2)/length(GhiaValues.v[:,1]))

            #Calculate maximum error per iteration
            pMaxError   = maximum(abs.(pOld.val - pCellVar.val))
            uMaxError   = maximum(abs.(uOld.val - uCellVar.val))
            vMaxError   = maximum(abs.(vOld.val - vCellVar.val))

            println("i=$i: u_RMSE=$uRMSE, v_RMSE=$vRMSE")
            println("Max Error: p=$pMaxError, u=$uMaxError, v=$vMaxError")
        end

         #Contour plot of lid driven cavity domain every 50 iterations
        if i%50 == 0
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

u, v = steadyLidDrivenCavity()

nx = u.domain.dims[1]
ny = v.domain.dims[2]
#Plot difference in reference data and numerical output
fig1 = plot(GhiaValues.u[:,1], GhiaValues.y*0.1,
            markershape = :circle,
            linetype    = :scatter,
            label       = "Ghia et al.")
       plot!(u.val[div(nx,2), 2:end-1], u.domain.cellCenters.y,
            label       = "Numerical Model",
            xlabel      = "x velocity (m/s)",
            ylabel      = "y (m)",
            legend      = :right)
display(fig1)

fig2 = plot(GhiaValues.x*0.1, GhiaValues.v[:,1],
            markershape = :circle,
            linetype    = :scatter,
            label       = "Ghia et al.")
       plot!(u.domain.cellCenters.y, v.val[2:end-1, div(ny,2)],
            label       = "Numerical Model",
            xlabel      = "x (m)",
            ylabel      = "y velocity (m)",
            legend      = :best)
display(fig2)
