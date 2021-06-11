function generateBC_1D(meshStruct::MeshStructure)
    #================================
    Generates a 1D boundary condition
    structure from a given mesh
    structure
    ================================#

    #Assume a 0 Neumann BC for left and right
    #of boundary

    left    = BC_Types()
    right   = BC_Types()

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
