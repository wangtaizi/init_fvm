function generateBC(meshStruct::MeshStructure)
    #==================================
    General function for formulation of
    Neuman boundary conditions by
    creating BC structure from
    mesh structure
    ==================================#

    dim = meshStruct.dimension

    if dim == 1
        BC = generateBC_1D(meshStruct)

    elseif dim == 2
        BC = generateBC_2D(meshStruct)
    end

    return BC
end

function generateBC_1D(meshStruct::MeshStructure)
    #================================
    Generates a 1D boundary condition
    structure from given mesh structure
    ================================#

    #Assume a 0 Neumann BC for left and right
    #of boundary

    left    = BC_Type()
    right   = BC_Type()

    left.neu = 1
    left.dir = 0
    left.val = 0

    right.neu = 1
    right.dir = 0
    right.val = 0

    #Generate empty arrays for other boundaries
    bottom  = [ ]
    top     = [ ]
    back    = [ ]
    front   = [ ]

    BC_1D = BoundaryCondition(meshStruct,
                left,
                right,
                bottom,
                top,
                back,
                front)

    return BC_1D
end

function generateBC_2D(meshStruct::MeshStructure)
    #================================
    Generates a 2D boundary condition
    structure from given mesh structure
    ================================#

    nx = meshStruct.dims[1]
    ny = meshStruct.dims[2]

    #Assume a 0 Neumann BC for left and right
    #of boundary
    left    = BC_Type()
    right   = BC_Type()
    top     = BC_Type()
    bottom  = BC_Type()

    left.neu = ones(1,ny)
    left.dir = zeros(1,ny)
    left.val = zeros(1,ny)

    right.neu = ones(1,ny)
    right.dir = zeros(1,ny)
    right.val = zeros(1,ny)

    top.neu = ones(nx,1)
    top.dir = zeros(nx,1)
    top.val = zeros(nx,1)

    bottom.neu = ones(nx,1)
    bottom.dir = zeros(nx,1)
    bottom.val = zeros(nx,1)

    #Generate empty arrays for other boundaries
    back    = [ ]
    front   = [ ]

    BC_2D = BoundaryCondition(meshStruct,
                left,
                right,
                bottom,
                top,
                back,
                front)

    return BC_2D
end
