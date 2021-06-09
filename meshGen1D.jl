module MeshGen1D

using ..MeshStructures

export meshGen1D

function meshGen1D(nx::Int, x_len::Float64)
    #===============================
    Generation of 1D uniform mesh
    nx: Total number of cells
    x_len: Length of numeric domain
    ===============================#

    dx = x_len / nx
    cellSize.x = dx*ones(nx+2, 1)
    cellSize.y = [0.0]
    cellSize.z = [0.0]

    cellCenter.x = (1:nx)*dx.-dx/2
    cellCenter.y = [0.0]
    cellCenter.z = [0.0]

    faceCenter.x = (0:nx)*dx
    faceCenter.y = [0.0]
    faceCenter.z = [0.0]

    mesh1D = MeshStructure(1,
                [nx, 1],
                cellSize,
                cellCenter,
                faceCenter,
                [1],
                [1])

    return mesh1D
end

function meshGen1D(xfaceCenters::Array)
    #==============================================
    Generation of 1D non-uniform mesh
    xfaceCenters: Array of x-dim cell face locations
    ==============================================#
    nx = length(xfaceCenters) - 1

    cellSize.x = [xfaceCenters[2]-xfaceCenters[1];
                    xfaceCenters[2:end]-xfaceCenters[1:end-1];
                    xfaceCenters[end]-xfaceCenters[end-1]]
    cellSize.y = [0.0]
    cellSize.z = [0.0]

    cellCenter.x = 0.5 * (xfaceCenters[2:end]+xfaceCenters[1:end-1])
    cellCenter.y = [0.0]
    cellCenter.z = [0.0]

    faceCenter.x = xfaceCenters
    faceCenter.y = [0.0]
    faceCenter.z = [0.0]

    mesh1D = MeshStructure(1,
                [nx, 1],
                cellSize,
                cellCenter,
                faceCenter,
                [1],
                [1])

    return mesh1D
end
