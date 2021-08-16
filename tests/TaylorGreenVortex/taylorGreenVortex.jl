using LinearAlgebra
using Plots

include("../../structs/structs.jl")
include("../../bcs/bcs.jl")
include("../../discret/discret.jl")
include("../../interpol/interpol.jl")
include("../../mesh/mesh.jl")
include("../../solvers/linearSolver.jl")
include("../../navier_stokes/calcCoeff.jl")
include("../../navier_stokes/momentum.jl")
include("../../navier_stokes/RhieChow.jl")

function taylorGreenVortex()
    #====================================================
    DESCRIPTION:
    The incompressible Navier-Stokes equations at the
    viscous limit are solved in a 2D boxed domain using
    the SIMPLE algorithm and Rhie-Chow interpolation
    method for a collocated grid.

    ∇⋅(μ∇ū) = ∇p
    ∇⋅ū = 0

    Numerical results are then compared with reference
    data given by Ghia et. al.:
     "High-Re solutions for incompressible flow using the Navier-Stokes
     equations and a multigrid method."
     Journal of computational physics 48.3 (1982): 387-411.

    =====================================================#
##
    #User Input Data
    pRelax      = 1.0                       #Pressure relaxation
    velRelax    = 0.9                       #Velocity relaxation
    uInit       = 0.0                       #Initial guess for x velocity
    vInit       = 0.0                       #Initial guess for y velocity
    pInit       = 0.0                       #Initial guess for pressure
    lidVel      = 1.0                       #Velocity of the lid
    Re          = 100                       #Reynolds number
    lx = ly     = pi                        #Length of square cavity
    rho         = 1000                      #Density of fluid
    mu          = rho*lx*lidVel/Re          #Dynamic viscosity
    ν           = mu/rho                    #Kinematic viscosity
    nIter       = 50                        #Number of iterations
    nx          = 30                        #Number of x nodes
    ny          = 30                        #Number of y nodes
    t_i         = 0.0
    t_f         = 0.5
    dt          = 0.05
    t           = t_i:dt:t_f
    x           = LinRange(0,lx,nx)
    y           = LinRange(0,ly,ny)
    F(t)        = exp(-2*ν*t)
    u(x,y,t)    = sin.(x).*cos.(y).*F.(t)
    v(x,y,t)    = -cos.(x).*sin.(y).*F.(t)


    #Formulate mesh
    msh = meshGen2D(nx, ny, lx, ly)

    #Create and set corresponding velocity boundary conditions
    uBC = generateBC(msh);  vBC = generateBC(msh)

    uBC.top.neu[:] .= 0;    uBC.top.dir[:] .= 1;    uBC.top.val[:] .= u(x,y[end],t_i)[:]
    uBC.bottom.neu[:] .= 0; uBC.bottom.dir[:] .= 1; uBC.bottom.val[:] .= u(x,y[1],t_i)[:]
    uBC.left.neu[:] .= 0;   uBC.left.dir[:] .= 1;   uBC.left.val[:] .= u(x[1],y,t_i)[:]
    uBC.right.neu[:] .= 0;  uBC.right.dir[:] .= 1;  uBC.right.val[:] .= u(x[end],y,t_i)[:]

    vBC.top.neu[:] .= 0;    vBC.top.dir[:] .= 1;    vBC.top.val[:] .= v(x,y[end],t_i)[:]
    vBC.bottom.neu[:] .= 0; vBC.bottom.dir[:] .= 1; vBC.bottom.val[:] .= v(x,y[1],t_i)[:]
    vBC.left.neu[:] .= 0;   vBC.left.dir[:] .= 1;   vBC.left.val[:] .= v(x[1],y,t_i)[:]
    vBC.right.neu[:] .= 0;  vBC.right.dir[:] .= 1;  vBC.right.val[:] .= v(x[end],y,t_i)[:]

    #Set the pressure correction equation boundary conditions
    pBC = generateBC(msh)           #Set to Neuman BC by default

    pBC.bottom.neu[div(nx,2)] = 0  #Set reference pressure value
    pBC.bottom.dir[div(nx,2)] = 1
    pBC.bottom.val[div(nx,2)] = 0

    #Build cell and face variables based off initial velocity and pressure guess
    muFaceVar   = generateFaceVar(msh, mu)
    rhoFaceVar  = generateFaceVar(msh, rho)
    faceVel     = generateFaceVar(msh, [uInit vInit])

    uCellVar    = generateCellVar(msh, uInit, uBC)
    vCellVar    = generateCellVar(msh, vInit, vBC)
    pCellVar    = generateCellVar(msh, pInit, pBC)

    #Take the gradient of the pressure
    pGradx, pGrady = grad(pCellVar)

    #main loop
    for j = 1:length(t)

        uBC.top.val[:] .= u(x,y[end],t[j])[:]
        uBC.bottom.val[:] .= u(x,y[1],t[j])[:]
        uBC.left.val[:] .= u(x[1],y,t[j])[:]
        uBC.right.val[:] .= u(x[end],y,t[j])[:]

        vBC.top.val[:] .= v(x,y[end],t[j])[:]
        vBC.bottom.val[:] .= v(x,y[1],t[j])[:]
        vBC.left.val[:] .= v(x[1],y,t[j])[:]
        vBC.right.val[:] .= v(x[end],y,t[j])[:]

        M_ut, RHS_ut = transient(uCellVar, dt)
        M_vt, RHS_vt = transient(vCellVar, dt)
        M_t          = M_ut + M_vt
        RHS_t        = RHS_ut + RHS_vt

        for i = 1:nIter
            #Define previous velocity iteration
            uOld = uCellVar
            vOld = vCellVar

            #Solve momentum equations
            uCellVar, vCellVar, u_ap, v_ap = momentum(msh, uOld, vOld, uBC, vBC,
                                                        pGradx, pGrady, faceVel,
                                                        rho, muFaceVar, velRelax,
                                                        "SIMPLE")

            #Pressure Correction Equation
            # ∇(Vₚ/aₚ)⋅∇p' = ∇⋅ū

            #Rhie Chow Interpolation to find Face Velocity
            faceVel = RhieChow(uCellVar, vCellVar, u_ap, v_ap, pCellVar)

            #Find the RHS value of the pressure correction by taking
            #divergence of the faceVel
            divFaceVel = divergence(faceVel)

            #face values of the velocity coefficients
            u_apFace = arithmeticFaceAvg(u_ap)
            v_apFace = arithmeticFaceAvg(v_ap)

            #Calculate the diffusion coefficient
            diffCoeffx = velRelax*msh.cellSize.y[1]./u_apFace.x
            diffCoeffy = velRelax*msh.cellSize.x[1]./v_apFace.y

            diffCoeff  = FaceVariable(msh, diffCoeffx, diffCoeffy, [])

            #Calculate diffusion term of pressure correction equation
            presCorrectionDiff = diffusionCD(diffCoeff)
            #diff = divergence(grad(diffCoeff))
            #println(presCorrectionDiff.x == diff)

            #Apply pressure boundary conditions
            M_pBC, RHS_pBC = applyBC(pBC)
            # There should be no need to impose pressure BCs to maintain continuity (at least in the fractional step method with a single pressure solve).
            # Instead the corrected face velocities should also satisfy the BCs.
            # Subtracting these from the original face velocities should give the correct boundary velocity corrections and these can be imposed.
            # The pressure field will then be floating with an additive constant, but otherwise the pressure gradients will be correct.

            #Solve for corrected pressure
            pNew = linearSolver(msh, presCorrectionDiff + M_pBC,
                                    divFaceVel[1] + RHS_pBC)

            #Update the old and new pressures and calculate max error
            pOld     = pCellVar
            pCellVar = CellVariable(pCellVar.domain, pCellVar.val + pNew.val*pRelax)

            #Update the velocity given new pressure
            pGradx, pGrady = grad(pNew)

            uNew = uCellVar.val - pGradx.domain.cellSize.y.*pGradx.val./u_ap.val

            vNew = vCellVar.val - pGrady.domain.cellSize.x.*pGrady.val./v_ap.val

            uCellVar = CellVariable(uCellVar.domain,
                        velRelax*uNew + (1-velRelax)*uCellVar.val)
            vCellVar = CellVariable(vCellVar.domain,
                        velRelax*vNew + (1-velRelax)*vCellVar.val)
            end



        #=
        #Contour plot of lid driven cavity domain every 50 iterations
        if j%div(length(t),3) == 0
            #pyplot()
            display(contourf(msh.cellCenters.x, msh.cellCenters.y,
                       hypot.(uCellVar.val[2:end-1,2:end-1]',
                       vCellVar.val[2:end-1,2:end-1]'),
                       ylims   = (0,pi),
                       xlims   = (0,pi),
                       #clims   = (0,1),
                       xlabel  = "x (m)",
                       ylabel  = "y (m)"))

        end
        =#

    end

    #Analytical solution at final timestep
    meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
    xx, yy        = meshgrid(x,y)
    uAnal   = u(xx,yy,t_f)
    vAnal   = v(xx,yy,t_f)

    return uCellVar, vCellVar, uAnal, vAnal, xx, yy
end

uCellVar, vCellVar, uAnal, vAnal, x, y  = taylorGreenVortex()

nx = uCellVar.domain.dims[1]
ny = vCellVar.domain.dims[2]

#Plot difference in analytical and numerical output
fig1 = plot(uCellVar.val[div(nx,2), 2:end-1], uCellVar.domain.cellCenters.y,
            markershape = :circle,
            linetype    = :scatter,
            label       = "Numerical Solution")
       plot!(reshape(uAnal,nx,ny)[div(nx,2), :], uCellVar.domain.cellCenters.y,
            label       = "Analytical Solution",
            xlabel      = "u (m/s)",
            ylabel      = "y (m)",
            legend      = :right)
display(fig1)

#Analytical vector field
fig2 = quiver(x,y, quiver=(0.1*uAnal, 0.1*vAnal))
display(fig2)

#Numerical vector field
fig3 = quiver(x,y, quiver=(0.1*uCellVar.val[:], 0.1*vCellVar.val[:]))
display(fig3)

#=
fig4 = plot(vCellVar.domain.cellCenters.y, vCellVar.val[2:end-1, div(ny,2)],
            markershape = :circle,
            linetype    = :scatter,
            label       = "Numerical Solution")
       plot!(vCellVar.domain.cellCenters.y, vAnal[:, div(ny,2)],
            label       = "Analytical Solution",
            xlabel      = "x (m)",
            ylabel      = "v (m)",
            legend      = :right)
display(fig4)
=#
