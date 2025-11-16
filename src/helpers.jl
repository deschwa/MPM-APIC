"""
(i  ,j+1)           (i+1,j+1)
        
           Particle

(i  ,j  )           (i+1,j  )
"""
function get_adjacent_grid_nodes(pos::AbstractVector{<:Real}, grid::Grid)

    i = floor(Int64, (pos[1]-grid.minx) / grid.dx) + 1  # 1-based indexing
    j = floor(Int64, (pos[2]-grid.miny) / grid.dy) + 1  # 1-based indexing

    return [
        (i, j),
        (i+1, j),
        (i, j+1),
        (i+1, j+1)
    ]
end



function cache_bin_map!(mp_group::MaterialPointGroup, grid::Grid)
    Nx = size(grid.mass, 1)
    Ny = size(grid.mass, 2)

    bin_map_cache = mp_group.bin_map_cache
    k_map_cache = mp_group.k_map_cache

    # Reset caches
    @inbounds for j in 1:Ny
        for i in 1:Nx
            empty!(bin_map_cache[i, j])
            empty!(k_map_cache[i, j])
        end
    end

    # Fill caches
    @inbounds for p_idx in 1:mp_group.N     
        for node_idx in 1:4
            i = mp_group.node_cache[1, node_idx, p_idx]
            j = mp_group.node_cache[2, node_idx, p_idx]
            push!(bin_map_cache[i, j], p_idx)
            push!(k_map_cache[i, j], node_idx)
        end
    end
end


function reset_grid!(grid::Grid)
    grid.v .= 0.0
    grid.v_new .= 0.0
    grid.momentum .= 0.0
    grid.momentum_new .= 0.0
    grid.f_ext .= 0.0
    grid.f_int .= 0.0
    grid.mass .= 0.0   
end