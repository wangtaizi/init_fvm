include("./structs/structs.jl")
include("./bcs/bcs.jl")
include("./discret/discret.jl")
include("./interpol/interpol.jl")
include("./mesh/mesh.jl")
include("./solvers/linearSolver.jl")

#Test 1D mesh domain
nx  = 4
lx  = 4
msh = meshGen1D(nx, lx)

#Test velocity and pressure field
u   = CellVariable(msh, [1.0 2.0 3.0 3.0])
v   = CellVariable(msh, zeros(1,4))
p   = CellVariable(msh, [0.6 0.8 0.7 0.6])

ap  = -1/16
uap = generateCellVar(msh, ap)
vap = CellVariable(msh, zeros(1,4))

rc  = RhieChow(u, v, uap, vap, p)
