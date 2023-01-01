
# module VoronoiTessellation

using Distances

# import GLMakie if you want to test or visualize your results
# using GLMakie
# using Plots

function make_neighbor_seeds_2D(p::Matrix{<:Number}, Lx::Number, Ly::Number)
    p1 = reshape(p[1,:],1,:)
    p2 = reshape(p[2,:],1,:)
    return [p [p1.+Lx; p2] [p1.+Lx; p2.+Ly] [p1; p2.+Ly] [p1.-Lx; p2] [p1.-Lx; p2.-Ly] [p1; p2.-Ly] [p1.+Lx; p2.-Ly] [p1.-Lx; p2.+Ly]]
end

function make_neighbor_seeds_3D(p::Matrix{<:Number}, Lx::Number, Ly::Number, Lz::Number)
    p1 = reshape(p[1,:],1,:)
    p2 = reshape(p[2,:],1,:)
    p3 = reshape(p[3,:],1,:)
    base = [p [p1.+Lx; p2; p3] [p1.+Lx; p2.+Ly; p3] [p1; p2.+Ly; p3] [p1.-Lx; p2; p3] [p1.-Lx; p2.-Ly; p3] [p1; p2.-Ly; p3] [p1.+Lx; p2.-Ly; p3] [p1.-Lx; p2.+Ly; p3]]
    base_add_Lz = copy(base)
    base_add_Lz[3,:] .+= Lz  
    base_minus_Lz = copy(base)
    base_minus_Lz[3,:] .-= Lz 
    return [base base_add_Lz base_minus_Lz]
end


function get_id_seeds(points::Matrix{<:Number}, seeds::Matrix{<:Number}, id_seeds::Vector{<:Int}; dist_func=Euclidean())
    
    num_points = size(points,2)
    range_points = 1:num_points
    num_seeds = length(id_seeds)
    range_seeds = 1:num_seeds
    id_seed_grids = ones(Int, num_points)

    for (id_p, p) in zip(range_points, eachcol(points))

        pair_dist = zeros(num_seeds)

        for (id_seed, seed) in zip(range_seeds, eachcol(seeds))
            @inbounds pair_dist[id_seed] = dist_func(p, seed)
        end

        (_, id_) = findmin(pair_dist)
        @inbounds id_seed_grids[id_p] = id_seeds[id_]

    end

    return id_seed_grids
end


function make_voro_3D(points::Matrix{<:Number}, num_seeds::Int; is_periodic=true, seeds=nothing, Lx=nothing, Ly=nothing, Lz=nothing)

    dim = 3
    if seeds===nothing
        seeds = rand(dim, num_seeds)
    end
    id_seeds = collect(1:num_seeds)

    # if periodic
    if is_periodic
        id_seeds_ = repeat(id_seeds,9*3)
        seeds_ = make_neighbor_seeds_3D(seeds, Lx+1e-3, Ly+1e-3, Lz+1e-3)
    else
        id_seeds_ = id_seeds
        seeds_ = seeds
    end

    mat_id = get_id_seeds(points, seeds_, id_seeds_)
    # mat_id_matrix = reshape(mat_id, nx, ny, nz);

    return mat_id

end


function test_voro_2D(nx::Int, ny::Int, Lx::Number, Ly::Number, num_seeds::Int; is_periodic=true, seeds=nothing)

    dim = 2
    if seeds===nothing
        seeds = rand(dim, num_seeds)
    end
    id_seeds = collect(1:num_seeds)

    # if periodic
    if is_periodic
        id_seeds_ = repeat(id_seeds, 9)
        seeds_ = make_neighbor_seeds_2D(seeds, Lx+1e-3, Ly+1e-3)
    else
        id_seeds_ = id_seeds
        seeds_ = seeds
    end

    x = collect(range(0,Lx,nx))
    y = collect(range(0,Ly,ny))
    X = [xx for xx in x, _ in y]
    Y = [yy for _ in x, yy in y]
    points = [X[:]'; Y[:]']

    mat_id = get_id_seeds(points, seeds_, id_seeds_)
    mat_id_matrix = reshape(mat_id, nx, ny);

    println("The Voronoi tessellation is done!")

    # heatmap(x, y, mat_id_matrix, aspect_ratio=:equal)
    # Plots.scatter!(seeds[1,:], seeds[2,:])
    # plt = Plots.scatter!(seeds_neighbors[1,:], seeds_neighbors[2,:])

    # To plot grain boundaries
    # Plots.contour(x, y, reshape(mat_id, 100, 100))

    fig = Figure()
    Axis(fig[1, 1], aspect=1.0)
    ms = meshscatter!(X[:], Y[:], color=mat_id_matrix[:] ,marker=Rect2((0,0), (1,1)), markersize=Lx/nx)
    Colorbar(fig[1,2], ms)
    resize_to_layout!(fig)

    return fig
end

function test_voro_3D(nx::Int, ny::Int, nz::Int, Lx::Number, Ly::Number, Lz::Number, num_seeds::Int; is_periodic=true, seeds=nothing)

    dim = 3
    if seeds===nothing
        seeds = rand(dim, num_seeds)
    end
    id_seeds = collect(1:num_seeds)

    # if periodic
    if is_periodic
        id_seeds_ = repeat(id_seeds,9*3)
        seeds_ = make_neighbor_seeds_3D(seeds, Lx+1e-3, Ly+1e-3, Lz+1e-3)
    else
        id_seeds_ = id_seeds
        seeds_ = seeds
    end

    x = collect(range(0,Lx,nx))
    y = collect(range(0,Ly,ny))
    z = collect(range(0,Lz,nz))
    X = [xx for xx in x, _ in y, _ in z]
    Y = [yy for _ in x, yy in y, _ in z]
    Z = [zz for _ in x, _ in y, zz in z]
    points = [X[:]'; Y[:]'; Z[:]']

    mat_id = get_id_seeds(points, seeds_, id_seeds_)
    mat_id_matrix = reshape(mat_id, nx, ny, nz);

    println("The Voronoi tessellation is done!")

    fig = Figure()
    Axis3(fig[1, 1], aspect=(1,1,1))
    ms = meshscatter!(X[:], Y[:], Z[:], color=mat_id_matrix[:] ,marker=Rect3((0,0,0), (1,1,1)), markersize=Lx/nx)
    Colorbar(fig[1,2], ms)
    resize_to_layout!(fig)
    return fig
end


# test_voro_2D(1000, 1000, 1.0, 1.0, 100, is_periodic=false)
# test_voro_3D(50, 50, 50, 1.0, 1.0, 1.0, 50, seeds=rand(3, 100))