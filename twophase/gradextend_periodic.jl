function gradextend_periodic(G::CellVariable)

	# extend result of grad from (nx,ny) to (nx+2,ny+2) by copying boundary values

	Gv = G.val
	nx = G.domain.dims[1]
	ny = G.domain.dims[2]
	Gv = Gv[2:nx+1,2:ny+1]

	#assemble matrix
	M1 = vcat(Gv[nx,ny], Gv[nx,:], Gv[nx,1])'
	M2 = hcat(Gv[:,ny], Gv, Gv[:,1])
	M3 = vcat(Gv[1,ny], Gv[1,:], Gv[1,1])'
	M = [M1; M2; M3]

	return(CellVariable(G.domain,M))	

end
