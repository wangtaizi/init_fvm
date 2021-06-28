function meshGen1D(nx, lx)
  #===============================
  Generation of 1D uniform mesh
  nx: Total number of cells
  lx: Length of numeric domain
  ===============================#
  cellSize    = Dim()
  cellCenter  = Dim()
  faceCenter  = Dim()

  dx = lx / nx

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
              (nx, 1),
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


function meshGen2D(nx, ny, lx, ly)
  #===============================
  Generation of 2D uniform mesh
  nx: Total number of x cells
  ny: Total number of y cells
  lx: Length of x domain
  ly: Length of y domain
  ===============================#
  cellSize    = Dim()
  cellCenter  = Dim()
  faceCenter  = Dim()

  dx = lx / nx
  dy = ly / ny

  cellSize.x = dx*ones(nx+2, 1)
  cellSize.y = dy*ones(ny+2, 1)
  cellSize.z = [0.0]

  cellCenter.x = (1:nx)*dx.-dx/2
  cellCenter.y = (1:ny)*dy.-dy/2
  cellCenter.z = [0.0]

  faceCenter.x = (0:nx)*dx
  faceCenter.y = (0:ny)*dy
  faceCenter.z = [0.0]

  cellNode     = reshape(1:(nx+2)*(ny+2), (nx+2, ny+2))
  cellNode     = cellNode[[1,end], [1,end]]

  mesh2D = MeshStructure(2,
              (nx, ny),
              cellSize,
              cellCenter,
              faceCenter,
              cellNode[:],
              [1])

  return mesh2D
end

function meshGen2D(xfaceCenters::Array, yfaceCenters::Array)
  #==============================================
  Generation of 2D non-uniform mesh
  xfaceCenters: Array of x-dim cell face locations
  yfaceCenters: Array of y-dim cell face locations
  ==============================================#


end
