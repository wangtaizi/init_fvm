using LinearAlgebra

include("./structs/structs.jl")
include("./bcs/bcs.jl")
include("./discret/discret.jl")
include("./interpol/interpol.jl")
include("./mesh/mesh.jl")
include("./solvers/linearSolver.jl")
include("./navier_stokes/calcCoeff.jl")
include("./navier_stokes/momentum.jl")
include("./navier_stokes/RhieChow.jl")

function steadyLidDrivenCavity()
    #====================================================
    DESCRIPTION:
    The incompressible Navier-Stokes equations are solved
    in a 2D boxed domain using the SIMPLE algorithm and
    Rhie-Chow interpolation method for a collocated grid.

    Numerical results are then compared with reference
    data given by Ghia et. al.:
     "High-Re solutions for incompressible flow using the Navier-Stokes
     equations and a multigrid method."
     Journal of computational physics 48.3 (1982): 387-411.

    =====================================================#

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
    nIter       = 250                       #Number of iterations
    nx          = 51                        #Number of x nodes
    ny          = 51                        #Number of y nodes

    #Formulate mesh
    msh = meshGen2D(nx, ny, cavityLen, cavityLen)

    #Create and set corresponding velocity boundary conditions
    uBC = generateBC(msh);  vBc = generateBC(msh)

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

    pBC.bottom.neu[end/2+0.5] .= 0  #Set reference pressure value
    pBC.bottom.dir[end/2+0.5] .= 1
    pBC.bottom.val[end/2+0.5] .= 0

    #Build cell and face variables based off initial velocity and pressure guess
    muFaceVar   = generateFaceVar(msh, mu)
    rhoFaceVar  = generateFaceVar(msh, rho)
    faceVel     = generateFaceVar(msh, [uInit vInit])

    uCellVar    = generateCellVar(msh, uInit, uBC)
    vCellVar    = generateCellVar(msh, vInit, vBC)
    pCellVar    = generateCellVar(msh, pInit, pBC)

    for 1:nIter
