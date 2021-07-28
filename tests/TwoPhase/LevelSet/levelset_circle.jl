include("../../../structs/structs.jl")
include("../../../mesh/mesh.jl")
include("../../../discret/discret.jl")
include("../../../twophase/twophase.jl")
include("../../../twophase/levelset.jl")


# create a 2D mesh

lx = 1
ly = 1
nx = 10
ny = 10
msh = meshGen2D(nx,ny,lx,ly)



# create a level set field of a circle

r = 0.2 #radius
x0 = 0.5 #x coordinate of circle center
y0 = 0.5 #y coordinate of circle center
dx = lx/nx
dy = ly/ny
x = repeat((0:nx+1)*dx.-dx/2,1,ny+2)
y = repeat((0:ny+1)'*dy.-dy/2,nx+2,1)
Gval = similar(x)
@. Gval = r - sqrt((x-x0)^2+(y-y0)^2)



# create level set, normal and curvature CellVariables

G = CellVariable(msh,Gval)
nx,ny = normal(G)
k = curvature(G)



# compute material variables
rho1 = 1000
rho2 = 1
mu1 = 1e-3
mu2 = 1e-5
F,rho,mu = levelset_materialproperties(G,rho1,rho2,mu1,mu2)
