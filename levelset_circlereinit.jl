include("./structs/structs.jl")
include("./mesh/mesh.jl")
include("./discret/discret.jl")
include("./twophase/twophase.jl")
include("./twophase/levelset.jl")


# create a 2D mesh

lx = 1
ly = 1
nx = 40
ny = 40
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

# declare material variables
rhoA = 1000
rhoB = 1
muA = 1e-3
muB = 1e-5

# create level set CellVariable

G = CellVariable(msh,Gval)



# compute quantities after 0 iterations

G0, Gs0 = levelset_reinit(G,0)
F0,rho0,mu0 = levelset_materialproperties(G0,rhoA,rhoB,muA,muB)
m0 = sum(F0.val)*dx*dy


# compute quantities after 1 iteration

G1, Gs1 = levelset_reinit(G0,1)
F1,rho1,mu1 = levelset_materialproperties(G1,rhoA,rhoB,muA,muB)
m1 = sum(F1.val)*dx*dy


# compute quantities after 2 iterations

G2, Gs2 = levelset_reinit(G1,1)
F2,rho2,mu2 = levelset_materialproperties(G2,rhoA,rhoB,muA,muB)
m2 = sum(F2.val)*dx*dy


# compute quantities after 2 iterations (at one go)

GX, GsX = levelset_reinit(G0,2)
FX,rhoX,muX = levelset_materialproperties(GX,rhoA,rhoB,muA,muB)
mX = sum(FX.val)*dx*dy
