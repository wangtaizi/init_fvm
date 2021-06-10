#===========================
Structure for elementary
Neumann boundary conditions
===========================#

struct BoundaryCondition
    domain::Any
    left::Any
    right::Any
    bottom::Any
    top::Any
    back::Any
    front::Any
end
