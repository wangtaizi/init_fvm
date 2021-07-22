function curvature(G::CellVariable)


	#Return CellVariable(s) corresponding to curvature associated with CellVariable G (first-order accuracy as curvature is computed at cell center and not interface)
	#G can be level set function or volume fraction field

	dim = G.domain.dimension

	
	if dim == 1

		kval = 0 #1D interface is flat by definition



	elseif dim == 2

		#Method 1: direct computation of surface Laplacian
		
		Gx, Gy = grad(G)
		Gx = gradextend(Gx)
		Gy = gradextend(Gy)
		Gxx, Gxy = grad(Gx)
		#Gxx = gradextend(Gxx)
		#Gxy = gradextend(Gxy)
		Gyx, Gyy = grad(Gy)
		#Gyx = gradextend(Gyx)
		#Gyy = gradextend(Gyy)
	
		Gxv = Gx.val
		Gyv = Gy.val
		Gxxv = Gxx.val
		Gxyv = Gxy.val	
		Gyyv = Gyy.val

		#use convention k = -div(n)
		kval = similar(Gxv)
		@. kval = (2. *Gxv*Gyv*Gxyv - Gyv*Gyv*Gxxv - Gxv*Gxv*Gyyv) / (Gxv*Gxv + Gyv*Gyv)^(3.0 / 2.0)

		
		#Method 2: computation from normal vector
		
		facenorm = faceGrad_2D(G)
		facenormx = facenorm.x
		facenormy = facenorm.y
		facenormn = similar(facenormx)
		@. facenormn = sqrt(facenormx.^2 + facenormy.^2)
		@. facenorm.x = facenorm.x / facenormn
		@. facenorm.y = facenorm.y / facenormn
		
		kvalface, divx, divy = divergence(facenorm)
		
		#Check that div grad = laplacian discretely
		sum(kvalface.val .- kval)
		
	end

	k = CellVariable(G.domain, kval)
	
	return k

end
