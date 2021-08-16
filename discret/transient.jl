function transient(oldϕ::CellVariable, dt::Float64, coeff::CellVariable)
    #====================================================
    DESCRIPTION:
    Discretizion of the transient term of the advective-
    diffusive equation using forward euler such that

    ∫ᵥ∂(ρϕ)/∂t dV ≃ (ρₚϕₚ - ρₚᵒˡᵈϕₚᵒˡᵈ)Vₚ / Δt

    where a constant coefficient (ρₚ) is assumed
    and provided by the user.

    RETURNS:
    Discretized matrix and RHS vector for the new
    time-stepped ϕ
    =====================================================#

    dim = oldϕ.domain.dimension

    if dim == 1
        M, RHS = transient1D(oldϕ::CellVariable, dt::Float64, coeff::CellVariable)
    elseif dim == 2
        M, RHS = transient2D(oldϕ::CellVariable, dt::Float64, coeff::CellVariable)
    end

    return M,RHS
end

function transient(oldϕ::CellVariable, dt::Float64)
    #====================================================
    DESCRIPTION:
    Discretizion of the transient term of the advective-
    diffusive equation using forward euler such that

    ∫ᵥ∂(ρϕ)/∂t dV ≃ (ρₚϕₚ - ρₚᵒˡᵈϕₚᵒˡᵈ)Vₚ / Δt

    where a constant coefficient (ρₚ) is assumed to be 1
    and the discretized equations are normalized by the
    cell volume Vₚ

    RETURNS:
    Discretized matrix and RHS vector for the new
    time-stepped ϕ
    =====================================================#

    dim = oldϕ.domain.dimension

    if dim == 1
        M, RHS = transient1D(oldϕ::CellVariable, dt::Float64)
    elseif dim == 2
        M, RHS = transient2D(oldϕ::CellVariable, dt::Float64)
    end

    return M,RHS
end

function transient1D(ϕold::CellVariable, dt::Float64)

    #mesh data
    nx    = ϕold.domain.dims[1]
    coeff = ones(nx,1)

    #build matrix data
    ind     = reshape(2:nx+1, nx, 1)
    diag    = reshape(coeff/dt, nx , 1)
    M       = sparse(ind, ind, diag, nx+2, nx+2)

    #build rhs data
    rhs      = zeros(nx+2, 1)
    rhs[ind] = reshape(coeff.*ϕold.val[2:nx+1]/dt, nx, 1)

    return M,rhs
end

function transient2D(ϕold::CellVariable, dt::Float64)

    #mesh data
    nx      = ϕold.domain.dims[1]
    ny      = ϕold.domain.dims[2]
    coeff   = ones(nx, ny)

    #build matrix data
    dom     = reshape(1:(nx+2)*(ny+2), nx+2, ny+2)
    ind     = reshape(dom[2:nx+1,2:ny+1], nx*ny, 1)
    diag    = reshape(coeff/dt, nx*ny , 1)
    M       = sparse(ind[:], ind[:], diag[:], (nx+2)*(ny+2), (nx+2)*(ny+2))

    #build rhs data
    rhs      = zeros((nx+2)*(ny+2), 1)
    rhs[ind] = reshape(coeff.*ϕold.val[2:nx+1,2:ny+1]/dt, nx*ny, 1)

    return M,rhs
end
