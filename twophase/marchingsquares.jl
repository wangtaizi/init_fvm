struct Cell
	crossings::Vector{UInt8}
end	

function marchingsquares(F::NodeVariable, G::CellVariable, XN::NodeVariable, YN::NodeVariable)

# Adapted from JuliaGeometry/Contour.jl @ GitHub

# The marching squares algorithm defines 16 cell types
# based on the edges that a contour line enters and exits
# through. The edges of the cells are identified using
# compass directions, while the vertices are ordered as
# follows:
#
#      N
#  4 +---+ 3
# W  |   |  E
#  1 +---+ 2
#      S
#
# Each cell type is identified with 4 bits, with each
# bit corresponding to a vertex (most significant bit -> 4, least significant bit -> 1).
# i.e. bit 4, bit 3, bit 2, bit 1
# A bit is set for vertex v_i is set if z(v_i) > 0. So a cell
# where a contour line only enters from the W edge and exits
# through the N edge will have the cell type: 0b0111
# Note that there are two cases where there are two
# lines crossing through the same cell: 0b0101, 0b1010.

N, S, E, W = (UInt8(1)), (UInt8(2)), (UInt8(4)), (UInt8(8))
NS, NE, NW = N|S, N|E, N|W
SN, SE, SW = S|N, S|E, S|W
EN, ES, EW = E|N, E|S, E|W
WN, WS, WE = W|N, W|S, W|E

# The way a contour crossing goes through a cell is labeled 
# by combining compass directions (e.g. a NW crossing connects
# the N and W edges of the cell). The Cell type records 
# the type of crossing that a cell contains. While most
# cells will have only one crossing, cell types 5 and 10 will
# have two crossings.

# Maps cell type to crossing types for non-ambiguous cells
edge_LUT = [SW, SE, EW, NE, 0, NS, NW, NW, NS, 0, NE, EW, SE, SW]

cells = Dict{(Tuple{Int,Int}),Cell}()

nx = F.domain.dims[1]
ny = F.domain.dims[2]

# Interested in zero contour
Fval = F.val
Gval = G.val
for i in 1:nx
	for j in 1:ny
		case = 1(Fval[i,j]>0) |
				2(Fval[i+1,j]>0) |
				4(Fval[i+1,j+1]>0) |
				8(Fval[i,j+1]>0)

		# contour does not go through these cells
		if case == 0 || case == 15
			continue
		end

		# process ambiguous cells (cases 5 and 10) using G as tiebreaker
		if case == 5
			if Gval[i,j] >= 0
				cells[(i,j)] = Cell([NW,SE])
			else
				cells[(i,j)] = Cell([NE,SW])
			end
		elseif case == 10
			if Gval[i,j] >= 0
				cells[(i,j)] = Cell([NE,SW])
			else
				cells[(i,j)] = Cell([NW,SE])
			end
		else
			cells[(i,j)] = Cell([edge_LUT[case]])
		end
	end
end	

# compute normal
NX, NY = normal(G)
nxval = NX.val
nyval = NY.val

# To be modified depending on desired output format: set up and output a list of points
# Second tuple: first integer tells us number of crossings, 
# followed by entries of unit normal vector,
# followed by coordinates of first (and second, if any) crossing. 
# Each crossing has two points (x1,y1,x2,y2), so four Float64s.
points = Dict{(Tuple{Int,Int}),(Tuple{Int,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64})}()
for (key,value) in cells
	crossing = first(value.crossings)
	edge1 = UInt8(0)
	edge2 = UInt8(0)
	for edge in [N,S,E,W]
		if edge & crossing != 0
			edge1 = edge
			edge2 = xor(crossing,edge1)
			break
		end
	end
	x1,y1 = interpolate_edge(XN.val,YN.val,Fval,key[1],key[2],edge1)
	x2,y2 = interpolate_edge(XN.val,YN.val,Fval,key[1],key[2],edge2)
	if length(value.crossings) == 1
		x3 = 0.0
		y3 = 0.0
		x4 = 0.0
		y4 = 0.0
	elseif length(value.crossings) == 2
		crossing = second(value.crossings)
		edge1 = UInt8(0)
		edge2 = UInt8(0)
		for edge in [N,S,E,W]
			if edge & crossing != 0
				edge1 = edge
				edge2 = xor(crossing,edge1)
				break
			end
		end
		x3,y3 = interpolate_edge(XN.val,YN.val,Fval,key[1],key[2],edge1)
		x4,y4 = interpolate_edge(XN.val,YN.val,Fval,key[1],key[2],edge2)
	end
	points[key] = (length(value.crossings),nxval[key[1]+1,key[2]+1],nyval[key[1]+1,key[2]+1],x1,y1,x2,y2,x3,y3,x4,y4)
end

return points

end

function interpolate_edge(x,y,z,xi::Int,yi::Int,edge::UInt8)
	N, S, E, W = (UInt8(1)), (UInt8(2)), (UInt8(4)), (UInt8(8))
	if edge == W
        y_interp = y[xi,yi] + (y[xi,yi + 1] - y[xi,yi]) * (-z[xi, yi]) / (z[xi, yi + 1] - z[xi, yi])
        x_interp = x[xi,yi]
    elseif edge == E
        y_interp = y[xi,yi] + (y[xi,yi + 1] - y[xi,yi]) * (-z[xi + 1, yi]) / (z[xi + 1, yi + 1] - z[xi + 1, yi])
        x_interp = x[xi + 1,yi]
    elseif edge == N
        y_interp = y[xi,yi + 1]
        x_interp = x[xi,yi] + (x[xi + 1,yi] - x[xi,yi]) * (-z[xi, yi + 1]) / (z[xi + 1, yi + 1] - z[xi, yi + 1])
    elseif edge == S
        y_interp = y[xi,yi]
        x_interp = x[xi,yi] + (x[xi + 1,yi] - x[xi,yi]) * (-z[xi, yi]) / (z[xi + 1, yi] - z[xi, yi])
    end
	
	return x_interp, y_interp
end
