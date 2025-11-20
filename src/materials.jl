using LinearAlgebra
using StaticArrays
using Base.Threads

# ---------------------------------
# Linear Elastic Isotropic Material
# ---------------------------------

function stress_update!(mp_group::MaterialPointGroup{LinearElastic}, dt::Float64)
    E = mp_group.material.E
    ν = mp_group.material.ν

    λ = (E * ν) / ((1 + ν) * (1 - 2 * ν))
    μ = E / (2 * (1 + ν))
    I_dim = @SMatrix [1.0 0.0;
                      0.0 1.0]
    
    @threads for p_idx in 1:mp_group.N
        Lp = @view mp_group.L[:,:,p_idx]
        σp = @view mp_group.σ[:,:,p_idx]

        ε_new = 0.5 * dt * (Lp + transpose(Lp))
        tr_ε = tr(ε_new)

        σ_new = σp + λ*tr_ε*I_dim + 2*μ*ε_new
        
        
        mp_group.σ[:,:,p_idx] .= σ_new
    end
end


