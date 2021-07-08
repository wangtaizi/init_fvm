include("./structs/structs.jl")
include("./bcs/bcs.jl")
include("./discret/discret.jl")
include("./interpol/interpol.jl")
include("./mesh/mesh.jl")
include("./solvers/linearSolver.jl")

#Test 1D mesh domain
nx  = 4
lx  = 1
ly  = 1
ny  = 1
msh = meshGen2D(nx, ny, lx, ly)

#Test velocity and pressure field
u   = CellVariable(msh, repeat([0 1.0 2.0 3.0 3.0 0]', 1, 3))
v   = CellVariable(msh, zeros(4,4))
p   = CellVariable(msh, repeat([0.6 0.8 0.7 0.6]', 1, 4))

ap  = -1/16
uap = generateCellVar(msh, ap)
vap = CellVariable(msh, zeros(4,4))

rc  = RhieChow(u, v, uap, vap, p)
