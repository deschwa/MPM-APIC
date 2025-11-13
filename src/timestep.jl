using StaticArrays
using Base.Threads



function timestep!(sim::MPMSimulation, alpha::Float64)
    grid = sim.grid
    Nx, Ny = size(grid.mass)
    mp_groups = sim.mp_groups

    
    # Reset grid
    reset_grid!(grid)

    # Cache Shape functions and bin maps for all particles
    @threads for mp_group in sim.mp_groups
        # Cache shape functions and nodes
        cache_shape_functions!(mp_group, grid)

        # cache bin maps
        cache_bin_map!(mp_group, grid)
    end

    # =========================
    # Particle to grid transfer
    # =========================
    @threads for i in 1:Nx
        for j in 1:Ny
            pos_ij = @SVector grid.pos[:,i,j]

            local_m = 0.0
            local_momentum = @MVector [0.0, 0.0]
            local_f_int = @MVector [0.0, 0.0]
            local_f_ext = @MVector [0.0, 0.0]

            for mp_group in mp_groups
                bin_map = mp_group.bin_map_cache
                
                for p_idx in bin_map[i,j]
                    mass = mp_group.mass[p_idx]

                    N_Ip = 0.0
                    ∇N_Ip = @SVector [0.0, 0.0]
                    for k in 1:4
                        if mp_group.node_cache[:,k,p_idx] .== (i,j)
                            N_Ip = mp_group.N_cache[k,p_idx]
                            ∇N_Ip = mp_group.∇N_cache[:,k,p_idx]
                            break
                        end # Node search
                    end #Node Search

                    # Transfer Mass
                    local_m += N_Ip * mass
                    
                    # Transfer momentum
                    local_momentum .+=  N_Ip * mass * (mp_group.vel[:,p_idx] 
                                        + mp_group.L[:,:,p_idx] * (mp_group.pos[:,p_idx] .- pos_ij)) # APIC Term

                    # calculate internal and external forces
                    local_f_ext .+= N_Ip * mp_group.ext_force_density[:,p_idx] * mass 
                    local_f_int .-= mp_group.volume[p_idx] * (mp_group.σ[:,:,p_idx] * ∇N_Ip)
                    

                end # loop over particles in corresponding cells               
            end # Loop over mp_groups

            # Set grid values
            grid.mass[i,j] = local_m
            grid.momentum[:,i,j] .= local_momentum
            grid.f_int[:,i,j] .= local_f_int
            grid.f_ext[:,i,j] .= local_f_ext

            # Update Momenta and velocities
            grid.momentum_new[:,i,j] .= grid.momentum[:,i,j] * (local_f_ext + local_f_int) * sim.dt
            if local_m > 1e-14
                grid.v[:,i,j] .= local_momentum / local_m
                grid.v_new[:,i,j] .= grid.momentum_new[:,i,j] / local_m
            end

        end # Loop over j ∈ 1:Ny
    end # Loop over i ∈ 1:Nx

    fix_boundaries!(grid)
    
    # =========================
    # Grid to Particle transfer
    # =========================
    for mp_group in mp_groups
        @threads for p_idx in 1:mp_group.N
            v_new_temp = @MVector [0.0, 0.0]
            x_new_temp = @MVector [0.0, 0.0]
            L_new_temp = @MMatrix [0.0 0.0; 0.0 0.0]
            for node_idx in 1:4
                i,j = mp_group.node_cache[:,node_idx,p_idx]
                N_Ip = mp_group.N_cache[node_idx, p_idx]
                ∇N_Ip = @SMatrix mp_group.∇N_cache[:, node_idx, p_idx]
                v_new_grid = @SVector grid.v_new[:,i,j]
                v_new_temp .+= N_Ip * (v_new_grid .- alpha.*grid.v[:,i,j])
                x_new_temp .+= N_Ip * v_new_grid
                L_new_temp .+= ∇N_Ip * v_new_grid'
            end

            mp_group.vel[:,p_idx] .*= alpha
            mp_group.vel[:,p_idx] .+= v_new_temp

            mp_group.pos[:,p_idx] .+= sim.dt * x_new_temp

            mp_group.L[:,:,p_idx] .= L_new_temp

            I_2 = @SMatrix [1.0 0.0; 0.0 1.0]
            mp_group.F[:,:,p_idx] .= (I_2 .+ alpha .* mp_group.L[:,:,p_idx]) * mp_group.F[:,:,p_idx]

            mp_group.volume[p_idx] = det(mp_group.F[:,:,p_idx]) * mp_group.volume_0[p_idx]
            

        end # Loop over particles

        stress_update!(mp_group, sim.dt)
    end # loop over mp_groups

    sim.t += sim.dt
end