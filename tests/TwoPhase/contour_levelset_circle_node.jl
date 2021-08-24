include("../../structs/structs.jl")
include("../../mesh/mesh.jl")
include("../../discret/discret.jl")
include("../../interpol/arithmeticNodeAvg.jl")
include("../../twophase/twophase.jl")


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



# create level set CellVariable and NodeVariable

G = CellVariable(msh,Gval)
F = arithmeticNodeAvg(G)

# create cell and node positions

XC = CellVariable(msh,x)
XN = arithmeticNodeAvg(XC)
YC = CellVariable(msh,y)
YN = arithmeticNodeAvg(YC)

# perform marching squares routine

cells = marchingsquares(F,G,XN,YN)
