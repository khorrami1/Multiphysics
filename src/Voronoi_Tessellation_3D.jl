
# module VoronoiTessellation

using Distances

func_dist(p1::Vector{Float64}, p2::Vector{Float64}) = sqrt((p1[1]-p2[1])^2+(p1[2]-p2[2])^2+(p1[3]-p2[3])^2)

# num_points = 100
dim = 3
# p = rand(num_points, dim)

num_seeds = 20
seeds = rand(num_seeds, dim)
# seeds = [0.24 0.5; 0.74 0.5; 0.5 0.74]
id_seeds = collect(1:num_seeds)

# id_seeds_neighbors = repeat(id_seeds,9)
id_seeds_neighbors = repeat(id_seeds,9*3) # for 3D

function make_neighbor_seeds_2D(p::Matrix{<:Number}, Lx::Number, Ly::Number)
    p1 = p[:,1]
    p2 = p[:,2]
    return [p ; p1.+Lx p2; p1.+Lx p2.+Ly ; p1 p2.+Ly ;p1.-Lx p2;
    p1.-Lx p2.-Ly ; p1 p2.-Ly; p1.+Lx p2.-Ly; p1.-Lx p2.+Ly]
end

function make_neighbor_seeds_3D(p::Matrix{<:Number}, Lx::Number, Ly::Number, Lz::Number)
    p1 = p[:,1]
    p2 = p[:,2]
    p3 = p[:,3]
    base = [p ; p1.+Lx p2 p3; p1.+Lx p2.+Ly p3; p1 p2.+Ly p3; p1.-Lx p2 p3;
    p1.-Lx p2.-Ly p3; p1 p2.-Ly p3; p1.+Lx p2.-Ly p3; p1.-Lx p2.+Ly p3]
    base_add_Lz = copy(base)
    base_add_Lz[:,3] .+= Lz  
    base_minus_Lz = copy(base)
    base_minus_Lz[:,3] .-= Lz 
    return [base; base_add_Lz; base_minus_Lz]
end

# func_periodic(x, L) = x - floor(x/L)*L

Lx = 1.0 + 0.0101
Ly = 1.0 + 0.0101
Lz = 1.0 + 0.0101

# seeds_neighbors = make_neighbor_seeds_2D(seeds, Lx, Ly)

# for 3D
seeds_neighbors = make_neighbor_seeds_3D(seeds, Lx, Ly, Lz)

func_pair_dist(p) = mapslices((p1)->Euclidean()(p1,p), seeds_neighbors, dims=2)

# custom_dist_func_2D(p1, p2) = sum((p1.-p2).^2)
# func_pair_dist(p) = mapslices((p1)->custom_dist_func_2D(p1,p), seeds_neighbors, dims=2)

# func_pair_dist(p[1,:])

function get_id_seeds(p::Vector{<:Number}, id_seeds_neighbors::Vector{<:Number})
    out = findmin(func_pair_dist(p))
    out2 = out[2]
    index = out2.I[1]
    return id_seeds_neighbors[index]
end


# mapslices(get_id_seeds, p, dims=1)

# ------------------------

nx = 50
ny = 50
nz = 50
x = collect(range(0,1,nx))
y = collect(range(0,1,ny))
z = collect(range(0,1,nz))

# for 2D
# X = [xx for xx in x, yy in y]
# Y = [yy for xx in x, yy in y]

# for 3D
X = [xx for xx in x, yy in y, zz in z]
Y = [yy for xx in x, yy in y, zz in z]
Z = [zz for xx in x, yy in y, zz in z]

# points = [X[:] Y[:]]

# for 3D
points = [X[:] Y[:] Z[:]]

mat_id = mapslices((p)->get_id_seeds(p, id_seeds_neighbors), points, dims=2)

using Plots

mat_id_matrix = reshape(mat_id, nx, ny, nz)

# heatmap(x, y, mat_id_matrix, aspect_ratio=:equal)
# Plots.scatter!(seeds[:,1], seeds[:,2])
# Plots.scatter!(seeds_neighbors[:,1], seeds_neighbors[:,2])

# To plot grain boundaries
# Plots.contour(x, y, reshape(mat_id, 100, 100))



using GLMakie

# data = ((i, j, k, RGBf(i, j, k)) for i in 0:0.1:1,
#                                      j in 0:0.1:1,
#                                      k in 0:0.1:1)

# x, y, z, color = (vec(getindex.(data, i)) for i in 1:4)

# meshscatter(x, y, z; color)


# fig = Figure(resolution=(800,800))
fig = Figure()
Axis3(fig[1, 1], aspect=(1,1,1))
ms = meshscatter!(X[:], Y[:], Z[:], color=mat_id_matrix[:] ,marker=Rect3((0,0,0), (1,1,1)), markersize=Lx/nx)
Colorbar(fig[1,2], ms)
resize_to_layout!(fig)

# end