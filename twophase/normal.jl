function normal(G::CellVariable)


	#Return CellVariable(s) corresponding to components of unit normal vector associated with CellVariable G
	#G can be level set function or volume fraction field

	dim = G.domain.dimension
	
	if dim == 1

		nx_orig = grad(G)
	
		nxval = similar(nx_orig.val)	
		@. nxval = nx_orig.val / abs(nx_orig.val) # should give +1 or -1 uniformly if G is monotonic

		nx = CellVariable(G.domain, nxval)

		return nx 


	elseif dim == 2

		nx_orig, ny_orig = grad(G)

		nval = similar(nx_orig.val)
		@. nval = sqrt(nx_orig.val^2 + ny_orig.val^2)

		nxval = similar(nx_orig.val)
		nyval = similar(ny_orig.val)
		@. nxval = nx_orig.val / nval
		@. nyval = ny_orig.val / nval

		nx = CellVariable(G.domain, nxval)
		ny = CellVariable(G.domain, nyval)
		
		return nx, ny

	
	end


end
