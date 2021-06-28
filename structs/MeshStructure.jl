#======================================
Structure of mesh and domain parameters
======================================#

mutable struct Dim
    x::Any
    y::Any
    z::Any

    Dim() = new()
end

struct MeshStructure
    dimension::Any
    dims::Any
    cellSize::Dim
    cellCenters::Dim
    faceCenters::Dim
    cornerNodes::Any
    cellEdges::Any

end
