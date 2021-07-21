function levelset_materialproperties(G::CellVariable, rho1, rho2, mu1, mu2)

	# Commented: Compute density (rho) and dynamic viscosity (mu) based on marker function (4.15) of Tryggvason, Scardovelli and Zaleski (2011) (c = 3)
	# Current version: Compute " " based on marker function (11) of van der Pijl, Segal, Vuik and Wesseling (2005) (c = 3/2)
	
	dim = G.domain.dimension

	if dim == 1

		a = similar(G.val)
		b = similar(G.val)
		@. a = abs(G.val) / G.domain.cellSize
		@. b = G.val / G.domain.cellSize
		

	elseif dim == 2

		nx = G.domain.dims[1]
		ny = G.domain.dims[2]
		dx = repeat(G.domain.cellSize.x,1,ny+2)
		dy = repeat(G.domain.cellSize.y',nx+2,1)
		ds = similar(dx)
		@. ds = sqrt(dx^2 + dy^2)

		a = similar(ds)
		b = similar(ds)	
		@. a = abs(G.val) / ds
		@. b = G.val / ds

	end

	#c = 3.
	c = 1.5

	# if a <= c then use interpolated marker function
	# if a > c and b > c then marker function = 1

	Fval = similar(G.val)
	rhoval = similar(G.val)
	muval = similar(G.val)
	#@. Fval = (a <= c) * (0.5 * (1 + b/c + sin(pi*b/c)/pi)) + (b > c)
	#@. rhoval = rho1 * Fval + rho2 * (1 - Fval)
	#@. muval = mu1 * Fval + mu2 * (1 - Fval)
	@. Fval = (a <= c) * (0.5 * (1 + sin(pi*b/2. / c))) + (b > c)
	@. rhoval = rho1 * Fval + rho2 * (1 - Fval)
	@. muval = mu1 * Fval + mu2 * (1 - Fval)
	
	F = CellVariable(G.domain,Fval)
	rho = CellVariable(G.domain,rhoval)
	mu = CellVariable(G.domain,muval)

	return F, rho, mu
		
end
