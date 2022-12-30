using Ferrite

# Generate a grid
N = 10
L = 1.0
dim = 3
left = zero(Vec{3})
right = L * ones(Vec{3})
grid = generate_grid(Tetrahedron, (N, N, N), left, right)

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