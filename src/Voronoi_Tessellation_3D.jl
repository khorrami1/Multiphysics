
using Distances

func_dist(p1::Vector{Float64}, p2::Vector{Float64}) = sqrt((p1[1]-p2[1])^2+(p1[2]-p2[2])^2+(p1[3]-p2[3])^2)

num_points = 100
p = rand(num_points, 2)

num_seeds = 3
seeds = rand(num_seeds, 2)
# seeds = [0.25 0.5; 0.75 1; 0.5 0.75]
id_seeds = collect(1:num_seeds)
id_seeds_neighbors = repeat(id_seeds,9)

function make_neighbor_seeds_2D(p::Matrix{<:Number}, Lx::Number, Ly::Number)
    p1 = p[:,1]
    p2 = p[:,2]
    return [p ; p1.+Lx p2; p1.+Lx p2.+Ly ; p1 p2.+Ly ;p1.-Lx p2;
    p1.-Lx p2.-Ly ; p1 p2.-Ly; p1.+Lx p2.-Ly; p1.-Lx p2.+Ly]
end

func_periodic(x, L) = x - floor(x/L)*L

Lx = 1.0 + 0.0101
Ly = 1.0 + 0.0101

seeds_neighbors = make_neighbor_seeds_2D(seeds, Lx, Ly)

func_pair_dist(p) = mapslices((p1)->SqEuclidean()(p1,p), seeds_neighbors, dims=2)
#func_pair_dist(p) = mapslices((p1)->custom_dist_func_2D(p1,p), seeds, dims=2)

func_pair_dist(p[1,:])

function get_id_seeds(p::Vector{<:Number}, id_seeds_neighbors::Vector{<:Number})
    out = findmin(func_pair_dist(p))
    out2 = out[2]
    index = out2.I[1]
    return id_seeds_neighbors[index]
end


# mapslices(get_id_seeds, p, dims=1)

# ------------------------
x = collect(range(0,1,100))
y = collect(range(0,1,100))

X = [xx for xx in x, yy in y]
Y = [yy for xx in x, yy in y]

points = [X[:] Y[:]]

mat_id = mapslices((p)->get_id_seeds(p, id_seeds_neighbors), points, dims=2)

using Plots

heatmap(x, y, reshape(mat_id, 100, 100))
# Plots.scatter!(seeds[:,1], seeds[:,2])
Plots.scatter!(seeds_neighbors[:,1], seeds_neighbors[:,2])

# To plot grain boundaries
Plots.contour(x, y, reshape(mat_id, 100, 100))



using GLMakie

data = ((i, j, k, RGBf(i, j, k)) for i in 0:0.1:1,
                                     j in 0:0.1:1,
                                     k in 0:0.1:1)

x, y, z, color = (vec(getindex.(data, i)) for i in 1:4)

meshscatter(x, y, z; color)



meshscatter(x, y, z; color, marker=Rect3((0,0,0), (1,1,1)))
