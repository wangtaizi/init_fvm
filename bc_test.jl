include("./structs/structs.jl")
include("./mesh/MeshGen.jl")
include("./bcs/generateBC.jl")
include("./bcs/applyBC.jl")

nx  = 10  #Number of cells
lx  = 1.0 #length of domain

msh = meshGen1D(nx, lx) #Generate 1D mesh
bc  = generateBC(msh)   #Generate bc structure
