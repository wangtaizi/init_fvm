module MeshStructures

#======================================
Structure of mesh and domain parameters
======================================#

export MeshStructure, cartesianCord

mutable struct cartesianCord
    x
    y
    z
    cartesianCord() = new()
end

struct MeshStructure
    dimension::Any
    dims::Any
    cellSize::cartesianCord
    cellCenters::cartesianCord
    faceCenters::cartesianCord
    cellNodes::Any
    cellEdges::Any

end

end
