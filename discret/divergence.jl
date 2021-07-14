function divergence(F::FaceVariable)
    #Given a field variable face vector,
    #return the divergence RHS vector

    dim = F.domain.dimension

    if dim == 1
        rhs_div = divergence_1D(F::FaceVariable)

        return rhs_div

    elseif dim == 2
        rhs_div, rhs_divx, rhs_divy = divergence_2D(F::FaceVariable)

        return rhs_div, rhs_divx, rhs_divy
    end
end

function divergence_1D(F::FaceVariable)
    #Calculate the 1D x gradient

end

function divergence_2D(F::FaceVariable)
    #calculate 2d gradient in x and y directions
    nx      = ϕ.domain.dims[1]
    ny      = ϕ.domain.dims[2]
    nodes   = reshape(1:(nx+2)*(ny+2), (nx+2,ny+2))
    xcells  = repeat(F.domain.cellSize.x[2:end-1], (1,ny))
    ycells  = repeat(F.domain.cellSize.y[2:end-1]', (nx,1))

    #Assign directional flux vectors
    Fe = F.x[2:nx+1,:]
    Fw = F.x[1:nx,:]
    Fn = F.y[:,2:ny+1]
    Fs = F.y[:,1+ny]

    divx = (Fe-Fw)./xcells
    divy = (Fn-Fs)./ycells

    #initialize RHS vector
    rhs_div  = zeros((nx+2)*(ny+2),1)
    rhs_divx = zeros((nx+2)*(ny+2),1)
    rhs_divy = zeros((nx+2)*(ny+2),1)

    #assign divergence to rhs vectors
    i = reshape(nodes[2:nx+1,2:ny+1], (nx*ny,1))

    rhs_div[i]   = reshape(divx+divy, (nx*ny,1))
    rhs_divx[i]  = reshape(divx, (nx*ny,1))
    rhs_divy[i]  = reshape(divy, (nx*ny,1))

    return rhs_div, rhs_divx, rhs_divy
end
