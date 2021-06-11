function diffusionCD(D::FaceVariable)
    #=======================================
    Given the diffusion coefficient, returns
    a discretized diffusion term using
    the method of central differences
    =======================================#

    dim = D.domain.dimension

    if d == 1
        M = diffusionCD_1D(D)
    end

    return M
end
