include("./mesh/MeshStructures.jl")
include("./mesh/MeshGen1D.jl")

using Plots
gr()

nx = 10;
L = 1;
ms = meshGen1D(nx, L)

f = plot(ms.cellCenters.x, ones(size(ms.cellCenters.x)),
        linecolor = :red, markershape = :circle,
        linetype = :scatter, markercolor = :transparent,
        label = "cell centers")
plot!(ms.faceCenters.x, ones(size(ms.faceCenters.x)),
        markershape = :cross,
        label =  "face centers",
        title = "1D Discretized Domain")
ylims!(0,2)
@show f
