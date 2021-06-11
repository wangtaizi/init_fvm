#===========================
Structure for elementary
boundary conditions and
their types
(Neumann, Dirichlet, Robin)
===========================#

mutable struct BC_Type
    neu
    dir
    val
    periodic::Bool

    BC_Type() = ( K = new(); K.periodic = false;
                return K)
end

struct BoundaryCondition
    domain::Any
    left::Any
    right::Any
    bottom::Any
    top::Any
    back::Any
    front::Any
end
