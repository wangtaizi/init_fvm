include("../../../structs/structs.jl")
include("../../../mesh/mesh.jl")
include("../../../bcs/bcs.jl")
include("../../../discret/discret.jl")
include("../../../interpol/arithmeticFaceAvg.jl")
include("../../../twophase/twophase.jl")
include("../../../twophase/levelset.jl")

#advection test for level set routine
case = 1

#material parameters
rho1 = 1000.	#density of phase 1
rho2 = 1.		#density of phase 2
mu1 = 1.0e-3	#dynamic viscosity of phase 1
mu2 = 1.0e-5	#dynamic viscosity of phase 2

#grid parameters
lx = 1
ly = 1
nx = 400
ny = 400

#create a 2D mesh and coordinates
msh = meshGen2D(nx,ny,lx,ly)
dx = lx/nx
dy = ly/ny
x = repeat((0:nx+1)*dx.-dx/2,1,ny+2)
y = repeat((0:ny+1)'*dy.-dy/2,nx+2,1)

#test cases
if case == 1 #circle translation

	#initialize interface
	r = 0.2 #radius
	x0 = 0.5 #x coordinate of circle center
	y0 = 0.5 #y coordinate of circle center
	Gval = similar(x)
	@. Gval = r - sqrt((x-x0)^2+(y-y0)^2)
	G = CellVariable(msh,Gval)

	#for verification
	G0 = G
	F0,rho0,mu0 = levelset_materialproperties(G0,rho1,rho2,mu1,mu2)

	#set velocity direction and magnitude
	Ux = 1. #x relative magnitude
	Uy = 2. #y relative magnitude
	Uval = 1. #velocity magnitude
	Ux = Ux/sqrt(Ux^2+Uy^2)*Uval
	Uy = Uy/sqrt(Ux^2+Uy^2)*Uval
	Ux = Ux * ones(nx+2,ny+2)
	Uy = Uy * ones(nx+2,ny+2)

	#advect

	#create periodic BC
	BC = generateBC(msh)
	BC.left.neu = zeros(1,ny)
	BC.right.neu = zeros(1,ny)
	BC.top.neu = zeros(nx,1)
	BC.bottom.neu = zeros(nx,1)
	BC.left.periodic = true
	BC.right.periodic = true
	BC.top.periodic = true
	BC.bottom.periodic = true

	#force G to be periodic
	ghost = ghostCells(Gval[2:nx+1,2:ny+1],BC)
	G = CellVariable(msh,ghost)
	G0 = G
	
	#create face variable for velocity and matrices for upwind convection
	Ucell = CellVariable(msh,Ux)
	Vcell = CellVariable(msh,Uy)
	Uface = arithmeticFaceAvg(Ucell)
	Vface = arithmeticFaceAvg(Vcell)
	U = FaceVariable(msh,Uface.x,Vface.y,[])
	M, Mx, My = upwindConvection(U)
	
	#set time step
	dt = 0.005/4

	#advect (to be written)
	for t = 1:1000*4
		Gval0 = G.val
		Gval0 = reshape(Gval0,((nx+2)*(ny+2),1))
		Gvalupdate = M * Gval0
		Gval0 = reshape(Gval0,(nx+2,ny+2))	
		Gvalupdate = reshape(Gvalupdate,(nx+2,ny+2))	
		Gval1 = similar(Gval0)
		@. Gval1 = Gval0 - dt * Gvalupdate
		global G = CellVariable(msh,Gval1)
	
		#force G to be periodic
		local ghost = ghostCells(Gval1[2:nx+1,2:ny+1],BC)
		global G = CellVariable(msh,ghost)
		
		if mod(t,4) == 0
			global G, Gs = levelset_reinit(G,BC,1)
		end
	end 
	

elseif case == 2 #Zalesak's disk (Zalesak (1999), see also Ansari (2019) Listing 3.1)

	#initialize interface
	Rval = similar(x)	
	@. Rval = sqrt((x-0.5)^2+(y-0.75)^2)
	Gval = similar(x)
	@. Gval = 0.15 - Rval
	bottom = 0.75 - 0.15 * cos(asin(0.025/0.15))
	for i = 1:nx+2
		for j = 1:ny+2
			Gnow = Gval[i,j]
			Rnow = Rval[i,j]
			xnow = x[i,j]
			ynow = y[i,j]
			if (ynow > 0.85 && Rnow < 0.15)
				Gnow = min(Gnow,ynow-0.85)
			end
			if (ynow < 0.85 && xnow < 0.475 && Rnow < 0.15)
				Gnow = min(Gnow,0.475-xnow)
			end
			if (ynow > 0.85 && xnow < 0.475 && Rnow < 0.15)
				Gnow = min(Gnow,sqrt((0.475-xnow)^2+(ynow-0.85)^2))
			end
			if (ynow < 0.85 && xnow > 0.525 && Rnow < 0.15)
				Gnow = min(Gnow,xnow-0.525)
			end
			if (ynow > 0.85 && xnow > 0.525 && Rnow < 0.15)
				Gnow = min(Gnow,sqrt((0.525-xnow)^2+(ynow-0.85)^2))
			end
			if (xnow > 0.475 && xnow < 0.525 && ynow < 0.85 && ynow > bottom)
				Gnow = min(0.85-ynow,min(0.525-xnow,xnow-0.475))
			end
			if (xnow > 0.475 && xnow < 0.525 && ynow < bottom)
				Gnow = min(sqrt((bottom-ynow)^2+(xnow-0.475)^2),sqrt((bottom-ynow)^2+(xnow-0.525)^2))
			end
			snow = Rnow < 0.15
			if (xnow > 0.475 && xnow < 0.525 && ynow < 0.85)
				snow = 0
			end
			snow = 2.0*snow - 1.0
			Gnow = abs(Gnow)*snow
			Gval[i,j] = Gnow
		end	
	end
	G = CellVariable(msh,Gval)
	BC = generateBC(msh) #default
	
	#for verification
	G0 = G
	F0,rho0,mu0 = levelset_materialproperties(G0,rho1,rho2,mu1,mu2)

	#set velocity field
	Ux = similar(x)
	Uy = similar(y)
	@. Ux = y - 0.5
	@. Uy = 0.5 - x

	#create face variable for velocity and matrices for upwind convection
	Ucell = CellVariable(msh,Ux)
	Vcell = CellVariable(msh,Uy)
	Uface = arithmeticFaceAvg(Ucell)
	Vface = arithmeticFaceAvg(Vcell)
	U = FaceVariable(msh,Uface.x,Vface.y,[])
	M, Mx, My = upwindConvection(U)

	#set time step
	dt = 2*pi/628/4 # equals dx (for 100x100 box) with some adjustment so it reaches exactly 2pi after 628 steps

	#advect (to be written)
	for t = 1:628*4
		Gval0 = G.val
		Gval0 = reshape(Gval0,((nx+2)*(ny+2),1))
		Gvalupdate = M * Gval0
		Gval0 = reshape(Gval0,(nx+2,ny+2))	
		Gvalupdate = reshape(Gvalupdate,(nx+2,ny+2))	
		Gval1 = similar(Gval0)
		@. Gval1 = Gval0 - dt * Gvalupdate
		global G = CellVariable(msh,Gval1)
		if mod(t,4) == 0
			global G, Gs = levelset_reinit(G,BC,1)
		end
	end 

end	

#for verification
F,rho,mu = levelset_materialproperties(G,rho1,rho2,mu1,mu2)

#plot

using Plots
default(show=false)
ENV["GKSwstype"]=100

xint = x[2:nx+1,1]
yint = y[1,2:ny+1]
G0int = G0.val[2:nx+1,2:ny+1]
p0 = Plots.contour(xint,yint,G0int',levels=[0.0],color=:black,colorbar=:none)
plotmain = Plots.plot(p0,aspect_ratio=:equal,xlim=[0,1],ylim=[0,1],xtickfont=font(20),ytickfont=font(20),grid=false,size=[1200,800])

if case == 1 
	savefig("interface_circle_init")
elseif case == 2
	savefig("interface_zalesak_init")
end

G1int = G.val[2:nx+1,2:ny+1]
p1 = Plots.contour!(xint,yint,G1int',levels=[0.0],color=:red,colorbar=:none)

if case == 1 
	savefig("interface_circle_fin")
elseif case == 2
	savefig("interface_zalesak_fin")
end

F0sum = sum(F0.val)*dx*dy
Fsum = sum(F.val)*dx*dy
print(abs(F0sum-Fsum)/F0sum*100)
