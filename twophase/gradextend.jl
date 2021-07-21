function gradextend(G::CellVariable)

	# extend result of grad from (nx,ny) to (nx+2,ny+2) by extrapolating with boundary gradients

	Gv = G.val
	nx = G.domain.dims[1]
	ny = G.domain.dims[2]

	#top of matrix
	N1 = Gv[1,:]
	N2 = Gv[2,:]
	N0 = 2*N1 - N2
	NW1 = N2[1]-N1[1]
	NE1 = N2[ny]-N1[ny]

	#bottom of matrix
	S1 = Gv[nx,:]
	S2 = Gv[nx-1,:]
	S0 = 2*S1 - S2
	SW1 = S2[1]-S1[1]
	SE1 = S2[ny]-S1[ny]

	#left of matrix
	W1 = Gv[:,1]
	W2 = Gv[:,2]
	W0 = 2*W1 - W2
	NW2 = W2[1]-W1[1]
	SW2 = W2[nx]-W1[nx]

	#right of matrix
	E1 = Gv[:,ny]
	E2 = Gv[:,ny-1]
	E0 = 2*E1 - E2
	NE2 = E2[1]-E1[1]
	SE2 = E2[nx]-E1[nx]

	#assemble matrix
	M1 = vcat(Gv[1,1]-NW1-NW2, N0, Gv[1,ny]-NE1-NE2)'
	M2 = hcat(W0, Gv, E0)
	M3 = vcat(Gv[nx,1]-SW1-SW2, S0, Gv[nx,ny]-SE1-SE2)'
	M = [M1; M2; M3]

	return(CellVariable(G.domain,M))	

end
