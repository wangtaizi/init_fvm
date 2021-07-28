include("../../../structs/structs.jl")
include("../../../mesh/mesh.jl")
include("../../../discret/discret.jl")
include("../../../twophase/twophase.jl")
include("../../../twophase/levelset.jl")

#advection test for level set routine
case = 2

#material parameters
rho1 = 1000.	#density of phase 1
rho2 = 1.		#density of phase 2
mu1 = 1.0e-3	#dynamic viscosity of phase 1
mu2 = 1.0e-5	#dynamic viscosity of phase 2

#grid parameters
lx = 1
ly = 1
nx = 40
ny = 40

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
	F,rho,mu = levelset_materialproperties(G,rho1,rho2,mu1,mu2)

	#set velocity direction and magnitude
	Ux = 1. #x relative magnitude
	Uy = 2. #y relative magnitude
	U = 1. #velocity magnitude
	Ux = Ux/sqrt(Ux^2+Uy^2)*U
	Uy = Uy/sqrt(Ux^2+Uy^2)*U

	#advect (to be written) (need periodic BC)



elseif case == 2 #Zalesak's disk (Zalesak (1999), see also Ansari (2019) Listing 3.1)

	#initialize interface
	Rval = similar(x)	
	@. Rval = sqrt((x-0.5)^2+(y-0.75)^2)
	Gval = similar(x)
	@. Gval = 0.15 - Rval
	bottom = 0.75 - 0.15 * cos(asin(0.025/0.15))
	for i = 1:nx+2
		for j = 1:nx+2
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
	
	#for verification
	F,rho,mu = levelset_materialproperties(G,rho1,rho2,mu1,mu2)

	#set velocity field
	Ux = similar(x)
	Uy = similar(y)
	@. Ux = y - 0.5
	@. Uy = 0.5 - x

	#set time step
	dt = 2*pi/628 # equals dx (for 100x100 box) with some adjustment so it reaches exactly 2pi after 628 steps

	#advect (to be written)


end	
