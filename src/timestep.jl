using StaticArrays
using Atomix: @atomic

const I_2 = @SMatrix [1.0 0.0; 0.0 1.0]




function timestep!(sim::MPMSimulation, alpha::Float64)
    grid = sim.grid
    Nx, Ny = size(grid.mass)
    mp_groups = sim.mp_groups

    
    # Reset grid
    reset_grid!(grid)

    # Cache Shape functions for all particles
    @threads for mp_group in mp_groups
        # Cache shape functions and nodes
        cache_shape_functions!(mp_group, grid)
    end


    # =========================
    # Particle to grid transfer
    # =========================
    for mp_group in mp_groups
        @threads for p_idx in 1:mp_group.N     

            m_p = mp_group.mass[p_idx]
            vol_p = mp_group.volume[p_idx]
            vel_p = SVector{2}(mp_group.vel[1,p_idx], mp_group.vel[2,p_idx])
            pos_p = SVector{2}(mp_group.pos[1,p_idx], mp_group.pos[2,p_idx])
            ext_acc_p = SVector{2}(mp_group.ext_acc[1,p_idx], mp_group.ext_acc[2,p_idx])
            L_p = SMatrix{2,2}(mp_group.L[1,1,p_idx], mp_group.L[2,1,p_idx],
                                mp_group.L[1,2,p_idx], mp_group.L[2,2,p_idx])
            σ_p = SMatrix{2,2}(mp_group.σ[1,1,p_idx], mp_group.σ[2,1,p_idx],
                                mp_group.σ[1,2,p_idx], mp_group.σ[2,2,p_idx])

            # affine_base = vel_p + L_p * pos_p

            for node_idx in 1:4
                i = mp_group.node_cache[1, node_idx, p_idx]
                j = mp_group.node_cache[2, node_idx, p_idx]

                if i >= 1 && i <= Nx && j >= 1 && j <= Ny
                    r_ij = SVector{2}(grid.pos[1,i,j], grid.pos[2,i,j])
    
                    N_Ip = mp_group.N_cache[node_idx, p_idx]
                    ∇N_Ip = SVector{2}(mp_group.∇N_cache[1, node_idx, p_idx], mp_group.∇N_cache[2, node_idx, p_idx])
    
                    # Transfer Mass
                    @atomic :monotonic grid.mass[i,j] += N_Ip * m_p
                    
                    # Transfer momentum
                    v_particle = vel_p - L_p * (pos_p - r_ij) # APIC Term
                    @atomic :monotonic grid.momentum[1,i,j] += N_Ip * m_p * v_particle[1]
                    @atomic :monotonic grid.momentum[2,i,j] += N_Ip * m_p * v_particle[2]
    
                    # calculate internal and external forces
                    @atomic :monotonic grid.f_ext[1,i,j] += N_Ip * ext_acc_p[1] * m_p 
                    @atomic :monotonic grid.f_ext[2,i,j] += N_Ip * ext_acc_p[2] * m_p 
    
                    f_int_particle = - vol_p * (σ_p * ∇N_Ip)
                    @atomic :monotonic grid.f_int[1,i,j] += f_int_particle[1]
                    @atomic :monotonic grid.f_int[2,i,j] += f_int_particle[2]
                end


            end # loop over nodes
        end # loop over particles
    end # loop over mp_groups
    
    
    # =================
    # Grid force update
    # =================
    @threads for i in 1:Nx
        for j in 1:Ny
            local_m = grid.mass[i,j]
            if local_m > 1e-14
                # Update Momenta and velocities
                grid.momentum_new[1,i,j] = grid.momentum[1,i,j] + (grid.f_ext[1,i,j] + grid.f_int[1,i,j]) * sim.dt
                grid.momentum_new[2,i,j] = grid.momentum[2,i,j] + (grid.f_ext[2,i,j] + grid.f_int[2,i,j]) * sim.dt
                
                grid.v[1,i,j] = grid.momentum[1,i,j] / local_m
                grid.v[2,i,j] = grid.momentum[2,i,j] / local_m
                
                grid.v_new[1,i,j] = grid.momentum_new[1,i,j] / local_m
                grid.v_new[2,i,j] = grid.momentum_new[2,i,j] / local_m
            end
        end
    end

    # ===================
    # Boundary conditions
    # ===================
    fix_boundaries!(grid)

    # =========================
    # Particle to grid transfer
    # =========================
    for mp_group in mp_groups
        @threads for p_idx in 1:mp_group.N
            v_new_temp = @MVector [0.0, 0.0]
            x_new_temp = @MVector [0.0, 0.0]
            L_new_temp = @MMatrix [0.0 0.0; 0.0 0.0]
            for node_idx in 1:4
                i,j = mp_group.node_cache[1,node_idx,p_idx], mp_group.node_cache[2,node_idx,p_idx]

                if i >= 1 && i <= Nx && j >= 1 && j <= Ny
                    N_Ip = mp_group.N_cache[node_idx, p_idx]
                    ∇N_Ip = SVector{2}(mp_group.∇N_cache[1, node_idx, p_idx], mp_group.∇N_cache[2, node_idx, p_idx])

                    v_grid = SVector{2}(grid.v[1,i,j], grid.v[2,i,j])
                    v_new_grid = SVector{2}(grid.v_new[1,i,j], grid.v_new[2,i,j])
                    v_new_temp .+= N_Ip * (v_new_grid - alpha*v_grid)

                    x_new_temp .+= N_Ip * v_new_grid
                    L_new_temp .+= v_new_grid * ∇N_Ip' 
                end
            end

            @views mp_group.vel[:,p_idx] .*= alpha
            @views mp_group.vel[:,p_idx] .+= v_new_temp

            @views mp_group.pos[:,p_idx] .+= sim.dt * x_new_temp

            @views mp_group.L[:,:,p_idx] .= L_new_temp * alpha

            F_p = SMatrix{2,2}(mp_group.F[1,1,p_idx], mp_group.F[2,1,p_idx],
                                mp_group.F[1,2,p_idx], mp_group.F[2,2,p_idx])
            L_p = SMatrix{2,2}(mp_group.L[1,1,p_idx], mp_group.L[2,1,p_idx],
                                mp_group.L[1,2,p_idx], mp_group.L[2,2,p_idx])

            F_new = (I_2 + sim.dt * L_p) * F_p
            @views mp_group.F[:,:,p_idx] .= F_new

            mp_group.volume[p_idx] = det(F_new) * mp_group.volume_0[p_idx]

        end # Loop over particles

        stress_update!(mp_group, sim.dt)
    end # loop over mp_groups

    sim.t += sim.dt
end