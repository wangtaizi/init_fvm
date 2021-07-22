function curvature(G::CellVariable)


	#Return CellVariable(s) corresponding to curvature associated with CellVariable G (first-order accuracy as curvature is computed at cell center and not interface)
	#G can be level set function or volume fraction field

	dim = G.domain.dimension

	
	if dim == 1

		kval = 0 #1D interface is flat by definition



	elseif dim == 2

		#ideally should use divergence function but divergence acts on a FaceVariable whereas components of the surface normal vector should be collocated
		#should verify that laplacian = div(grad) discretely but perhaps these are not the right variables to use for verification

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

	end

	k = CellVariable(G.domain, kval)
	
	return k

end
