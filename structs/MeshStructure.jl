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
    dimension::Any      #Number of dimension of domain
    dims::Any           #Domain shape
    cellSize::Dim       #Size of each cell
    cellCenters::Dim    #Location of cell centers
    faceCenters::Dim    #location of face centers
    cornerNodes::Any    #Vector of corner nodes
    cellEdges::Any      #Cell edges 

end
