function levelset_reinit(G::CellVariable, BC::BoundaryCondition, n)

	#reinitialize G to make sure it remains as close to a distance function as possible
	#performed as necessary (every time step, or every few time steps) with n iterations
	#reinitialization shifts the G=0 isosurface and causes mass loss - thus should be sparingly used 
	#see Gomez, Hernandez, Lopez (2005) and Herrmann (2008) for thresholds on |grad G| to trigger this routine

	#calculate cell size
	nx = G.domain.dims[1]
	ny = G.domain.dims[2]
	dx = repeat(G.domain.cellSize.x,1,ny+2)
	dy = repeat(G.domain.cellSize.y',nx+2,1)
	ds = similar(dx)
	@. ds = sqrt(dx^2 + dy^2)

	#set up sign function as per eq (36) of Peng, Merriman, Osher, Zhao and Kang (1999)
	Gv = G.val
	Gx, Gy = grad(G)
	if (BC.left.periodic == true && BC.right.periodic == true && BC.top.periodic == true && BC.bottom.periodic == true)
		Gx = gradextend_periodic(Gx)
		Gy = gradextend_periodic(Gy)
	else
		Gx = gradextend(Gx)
		Gy = gradextend(Gy)
	end
	Gxv = Gx.val
	Gyv = Gy.val
	Gsv = similar(Gxv)
	@. Gsv = sqrt(Gxv^2 + Gyv^2) #magnitude of grad G
	S = similar(Gsv)
	@. S = Gv / sqrt(Gv^2 + Gsv^2 * ds^2)

	#the Eikonal equation need only be solved near the interface using a fast marching method
	#but we can do it everywhere for now with a Hamilton-Jacobi PDE-based scheme as long as the mesh is reasonably small

	#Higher order TVD schemes may be used with high-order TVD RK time advancement as per Peng et al. (1999)
	#But we use the current second-order grad operators (as above) with forward Euler for now for simplicity

	for i = 1:n

		#select time step using (39) of Peng et al. (1999)
		dt = 0.49 * minimum(ds) / maximum(abs.(S))

		# Solve eq (37) of Peng et al. (1999)
		Gvnew = similar(Gv)
		@. Gvnew = Gv - dt * S * (Gsv - 1.)

		# update variables
		Gv = Gvnew
		Gx, Gy = grad(CellVariable(G.domain,Gv))
		if (BC.left.periodic == true && BC.right.periodic == true && BC.top.periodic == true && BC.bottom.periodic == true)
			Gx = gradextend_periodic(Gx)
			Gy = gradextend_periodic(Gy)
		else
			Gx = gradextend(Gx)
			Gy = gradextend(Gy)
		end
		Gxv = Gx.val
		Gyv = Gy.val
		Gsv = similar(Gxv)
		@. Gsv = sqrt(Gxv^2 + Gyv^2) #magnitude of grad G
		S = similar(Gsv)
		@. S = Gv / sqrt(Gv^2 + Gsv^2 * ds^2)

	end	

	Gnew = CellVariable(G.domain,Gv)
	Gsnew = CellVariable(G.domain,Gsv) #used for debugging, can comment out later
	
	return Gnew, Gsnew

end
