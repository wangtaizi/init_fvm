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

function poiseuilleFlow()
    #====================================================
    DESCRIPTION:
    The incompressible Navier-Stokes equations are solved
    for a fully developed cacnonical poiseuille flow
    within a 2D channel
    =====================================================#
##
    #User Input Data
    velRelax    = 1.0                       #Velocity relaxation
    uInit       = 0.0                       #Initial guess for x velocity
    vInit       = 0.0                       #Initial guess for y velocity
    pInit       = 0.0                       #Initial guess for pressure
    height      = 1.0                       #channel height
    len         = 2*height                  #channel length
    dPdx        = -0.01                      #x pressure gradient
    rho         = 1000                      #Density of fluid
    mu          = .001                      #Dynamic viscosity
    nx          = 500                       #Number of x nodes
    ny          = 500                       #Number of y nodes

    #Formulate mesh
    msh = meshGen2D(nx, ny, len, height)

    #Create and set corresponding velocity boundary conditions
    uBC = generateBC(msh);  vBC = generateBC(msh)

    uBC.top.neu[:]    .= 0; uBC.top.dir[:]    .= 1; uBC.top.val[:]    .= 0 #surface velocity
    uBC.bottom.neu[:] .= 0; uBC.bottom.dir[:] .= 1; uBC.bottom.val[:] .= 0 #No slip
    uBC.left.neu[:]   .= 1; uBC.left.dir[:]   .= 0; uBC.left.val[:]   .= 0 #No penetration
    uBC.right.neu[:]  .= 1; uBC.right.dir[:]  .= 0; uBC.right.val[:]  .= 0 #No penetration

    vBC.top.neu[:]    .= 0; vBC.top.dir[:]    .= 1; vBC.top.val[:]    .= 0 #No penetration
    vBC.bottom.neu[:] .= 0; vBC.bottom.dir[:] .= 1; vBC.bottom.val[:] .= 0 #No penetration
    vBC.left.neu[:]   .= 0; vBC.left.dir[:]   .= 1; vBC.left.val[:]   .= 0 #No slip
    vBC.right.neu[:]  .= 0; vBC.right.dir[:]  .= 1; vBC.right.val[:]  .= 0 #No slip

    #Set the pressure correction equation boundary conditions
    pBC = generateBC(msh)           #Set to Neuman BC by default

    pBC.top.neu[:]    .= 1; pBC.top.dir[:]    .= 0; pBC.top.val[:]    .= 0
    pBC.bottom.neu[:] .= 1; pBC.bottom.dir[:] .= 0; pBC.bottom.val[:] .= 0
    pBC.left.neu[:]   .= 1; pBC.left.dir[:]   .= 0; pBC.left.val[:]   .= 0
    pBC.right.neu[:]  .= 1; pBC.right.dir[:]  .= 0; pBC.right.val[:]  .= 0

    pBC.bottom.neu[div(ny,2)] = 0; pBC.bottom.dir[div(ny,2)] = 1; pBC.bottom.val[div(ny,2)] = 1

    #Build cell and face variables based off initial velocity and pressure guess
    muFaceVar   = generateFaceVar(msh, mu)
    faceVel     = generateFaceVar(msh, [uInit vInit])

    uCellVar    = generateCellVar(msh, uInit, uBC)
    vCellVar    = generateCellVar(msh, vInit, vBC)
    pCellVar    = generateCellVar(msh, 1, pBC)

    pGradx      = generateCellVar(msh, dPdx, pBC)
    pGrady      = generateCellVar(msh, 0.0, pBC)

    #Solve momentum equations
    uCellVar, vCellVar, u_ap, v_ap = momentum(msh, uCellVar, vCellVar, uBC, vBC,
                                                pGradx, pGrady, faceVel, rho,
                                                muFaceVar, velRelax, "SIMPLE")


    #Analytical solution
    U(y) = (-dPdx/(2*mu)).*y.*(height.-y)
    yRng = LinRange(0, height, ny)
    uAnalytical = U(yRng)

    #Calculate root mean squared error
    uRMSE   = sqrt(sum((uAnalytical[:].-uCellVar.val[3,2:end-1]).^2)/length(uAnalytical[:]))
    println(" ")
    println("u_RMSE=$uRMSE")

    #Plot difference in reference data and numerical output
    fig1 = plot(uAnalytical[:], yRng,
                #markershape = :circle,
                linestyle   = :dash,
                linewidth   = 3,
                label       = "Analytical Solution")
           plot!(uCellVar.val[3, 2:end-1], uCellVar.domain.cellCenters.y,
                label       = "Numerical Solution",
                xlabel      = "u (m/s)",
                ylabel      = "y (m)",
                legend      = :topright)
    display(fig1)
end

poiseuilleFlow()
