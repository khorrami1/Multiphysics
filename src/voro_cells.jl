using Ferrite
include("Voronoi_Tessellation.jl")

# Generate a grid
N = 50
L = 1.0
dim = 3
left = zero(Vec{3})
right = L * ones(Vec{3})
cell_type = QuadraticTetrahedron 
grid = generate_grid(Cell{3,20,6}, (N, N, N), left, right)

cells = grid.cells
nodes = grid.nodes

num_cells = getncells(grid)
centers_cell = zeros(dim, num_cells)

for id_cell in 1:num_cells
    coords_cell = getcoordinates(grid, id_cell)
    centers_cell[:, id_cell] = sum(coords_cell)/length(coords_cell)
end

centers_cell

# using FerriteViz
# import GLMakie
# FerriteViz.wireframe(grid,markersize=10,strokewidth=2)

# mat_id = make_voro_3D(centers_cell, 5, is_periodic=false)

mat_id = make_voro_3D(centers_cell, 50, is_periodic=true, Lx=L+1e-3, Ly=L+1e-3, Lz=L+1e-3)

# vtk_file = vtk_grid("voro_3D", grid)
# vtk_cell_data(vtk_file, mat_id, "mat_id")
# close(vtk_file)

vtk_grid("voro_3D", grid) do vtk_file
    vtk_cell_data(vtk_file, mat_id, "mat_id")
end