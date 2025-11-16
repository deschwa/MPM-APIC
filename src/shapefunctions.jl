using Base.Threads
using LinearAlgebra
using StaticArrays



function shape_function(r_rel::SVector{2,Float64}, dx::Float64, dy::Float64)
           
    N_x = max(1 - abs(r_rel[1]) / dx, 0.0)
    N_y = max(1 - abs(r_rel[2]) / dy, 0.0)
    N_I = N_x * N_y


    dN_xdx = - sign(r_rel[1]) / dx
    dN_ydy = - sign(r_rel[2]) / dy
    ∇N_I = @SVector [dN_xdx * N_y,
            dN_ydy * N_x]


    return (N_I, ∇N_I)
end


function cache_shape_functions!(mp_group::MaterialPointGroup, grid::Grid)

    @threads for p_idx in 1:mp_group.N
        @views pos_p = mp_group.pos[:, p_idx]

        adj_nodes = get_adjacent_grid_nodes(pos_p, grid)

        for (node_idx, (i,j)) in enumerate(adj_nodes)
            mp_group.node_cache[:, node_idx, p_idx] .= (i,j)

            rel_pos = SVector{2}(pos_p .- grid.pos[:, i, j])

            N_Ip, ∇N_Ip = shape_function(rel_pos, grid.dx, grid.dy)

            mp_group.N_cache[node_idx, p_idx] = N_Ip
            mp_group.∇N_cache[:, node_idx, p_idx] .= ∇N_Ip

        end

    end
    
end