#===========================
Structure for elementary
boundary conditions and
their types
(Neumann, Dirichlet, Robin)
===========================#

mutable struct BC_Type
    neu #Neumann Condition
    dir #Dirichlet Condition
    val #Value of Condition
    periodic::Bool #Used to check whether BC is periodic

    BC_Type() = ( K = new(); K.periodic = false;
                return K)
end

struct BoundaryCondition
    domain::MeshStructure
    left::Any
    right::Any
    bottom::Any
    top::Any
    back::Any
    front::Any
end
