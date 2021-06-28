include("./structs/structs.jl")
include("./mesh/meshGen.jl")
include("./mesh/generateCellVar.jl")
include("./bcs/generateBC.jl")
include("./bcs/applyBC.jl")
include("./bcs/ghostCells.jl")
include("./interpol/cellFaceHarmonic.jl")
include("./discret/diffusionCD.jl")
include("./solvers/linearSolver.jl")

#Solving Diffusion equation
#∇⋅(-D∇c) = 0

lx  = 0.01               #domain  length
nx  = 10                #cell total number
msh = meshGen1D(nx, lx) #mesh gen
bc  = generateBC(msh)   #bc structure gen

#Applying dirichlet BCs at bounds
bc.left.neu     = 0 #No neuman condition
bc.left.dir     = 1 #Use dirichlet condition
bc.left.val     = 1 #c = 1 @ x = 0

bc.right.neu    = 0 #No neuman
bc.right.dir    = 1 #Use dirichlet
bc.right.val    = 0 # c = 0 @ x = L

#Define diffusivity coefficient
D       = generateCellVar(msh, 1e-5) #assign const val to each cell
D_face  = cellFaceHarmonic(D)        #assign harmonic mean face vals

#Formulate coefficient matrix and solve
M_diff          = diffusionCD(D_face)
M_bc, RHS_bc    = applyBC(bc)

sol             = linearSolver(msh, M_diff + M_bc, RHS_bc)

#Plot 1D solution
println(length(sol.val))
using Plots

f = plot(LinRange(0, lx, length(sol.val)), sol.val)
@show f
