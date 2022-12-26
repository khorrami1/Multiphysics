
func_dist(p1::Vector{Float64}, p2::Vector{Float64}) = sqrt((p1[1]-p2[1])^2+(p1[2]-p2[2])^2+(p1[3]-p2[3])^2)

num_points = 100
p = rand(3, num_points)

num_seeds = 3
seeds = rand(3, num_seeds)

func_pair_dist(p) = mapslices((p1)->func_dist(p1,p), seeds, dims=1)

func_pair_dist(p[:,1])

function get_id_seeds(p::Vector{Float64})
    out = findmin(func_pair_dist(p))
    out2 = out[2]
    return out2.I[2]
end


mapslices(get_id_seeds, p, dims=1)