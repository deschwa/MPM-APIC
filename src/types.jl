using StaticArrays



"""
Materials
"""
abstract type AbstractMaterial end

struct LinearElastic<:AbstractMaterial
    E::Float64        # Young's Modulus
    ν::Float64        # Poisson's Ratio
    ρ::Float64        # Density
end

struct NeoHookean<:AbstractMaterial
    μ::Float64        # Shear Modulus
    λ::Float64       # Lamé's First Parameter
    ρ::Float64        # Density
end




"""
Material Points
"""
struct MaterialPointGroup{MaterialType<:AbstractMaterial}
    N::Int64                        # Nr of Particles in Group
    
    # Vector quantities. Access Vectors via pos_p = pos[:, particle_idx] and similar
    pos::Array{Float64, 2}          # Positions
    vel::Array{Float64, 2}          # Velocities
    ext_acc::Array{Float64, 2}      # External Forces

    # Scalar quantities
    mass::Vector{Float64}           # Masses
    volume::Vector{Float64}         # Volumes
    volume_0::Vector{Float64}       # Initial Volumes
    density::Vector{Float64}        # Densities

    # Matrix quantities. Access Matrices via F_p = F[:, :, particle_idx] and similar
    F::Array{Float64, 3}       # Deformation Gradients
    σ::Array{Float64, 3}       # Cauchy Stresses  
    L::Array{Float64, 3}       # Velocity Gradients

    material::MaterialType                  # Material Model
    type::String                            # Material Type

    node_cache::Array{Int64, 3}                 # node cache; 2x4xN Matrix, access via i, j = node_cache[:, grid_index, particle_idx]
    N_cache::Array{Float64, 2}                  # N_Ip cache; 4xN Matrix, access via N_Ip = N[node_idx, particle_idx]
    ∇N_cache::Array{Float64, 3}                 # ∇N Cache; 2x4xN Matrix, access via ∇N_Ip = ∇N_cache[:, node_idx, particle_idx]
    
    # Constructor using only pos, vel, mass, volume, and material. AbstractString is used because sometimes CSVs are read as String7
    function MaterialPointGroup(pos::Array{Float64, 2}, 
                                vel::Array{Float64, 2},
                                mass::Vector{Float64}, 
                                volume::Vector{Float64}, 
                                material::MaterialType, type::AbstractString) where {MaterialType<:AbstractMaterial}

        # Number of particles is the number of columns in pos (pos is 2 x N)
        N = size(pos, 2)

        # Basic consistency checks to avoid bounds errors later
        @assert size(pos, 1) == 2 "pos must be a 2 x N array"
        @assert size(vel, 1) == 2 "vel must be a 2 x N array"
        @assert size(vel, 2) == N "vel must have the same number of columns as pos"
        @assert length(mass) == N "mass length must equal number of particles"
        @assert length(volume) == N "volume length must equal number of particles"

        density = [m/v for (m,v) in zip(mass, volume)]
        volume_0 = copy(volume)

        F = zeros(Float64, (2,2,N))
        for p_idx in 1:N
            F[:,:,p_idx] .= [1.0 0.0; 0.0 1.0]
        end

        σ = zeros(Float64, (2,2,N))
        L = zeros(Float64, (2,2,N))

        ext_acc = zeros(Float64, (2,N))

        node_cache = zeros(Int64, (2,4,N))
        N_cache = zeros(Float64, (4,N))
        ∇N_cache = zeros(Float64, (2,4,N))

        # Ensure stored `type` is a concrete String to match the field type
        type_str = String(type)

        new{MaterialType}(N,
            pos,
            vel,
            ext_acc,
            mass,
            volume,
            volume_0,
            density,
            F,
            σ,
            L,
            material,
            type_str,
            node_cache,
            N_cache,
            ∇N_cache)

    end

end


"""
struct Grid

Contains all necessary information about the grid.
"""
struct Grid
    pos::Array{Float64, 3}                  # (2)x(Nx)x(Ny) Array with Positions

    v::Array{Float64, 3}                    # (2)x(Nx)x(Ny) Array with Velocities
    v_new::Array{Float64, 3}                # (2)x(Nx)x(Ny) Array with New Velocities

    momentum::Array{Float64, 3}             # (2)x(Nx)x(Ny) Array with Momenta
    momentum_new::Array{Float64, 3}         # (2)x(Nx)x(Ny) Array with New Momenta

    f_ext::Array{Float64, 3}                # (2)x(Nx)x(Ny) Array with External Forces
    f_int::Array{Float64, 3}                # (2)x(Nx)x(Ny) Array with Internal Forces

    mass::Array{Float64, 2}                 # (Nx)x(Ny) Array with Masses

    dx::Float64                             # Grid Spacing in x
    dy::Float64                             # Grid Spacing in y
    minx::Float64                           # Minimum x-coordinate
    miny::Float64                           # Minimum y-coordinate

    function Grid(Nx::Int64, Ny::Int64, minx::Float64, maxx::Float64, miny::Float64, maxy::Float64)
        dx::Float64 = (maxx - minx) / (Nx-1)
        dy::Float64 = (maxy - miny) / (Ny-1)

        xs = round.(collect(range(minx, stop=maxx, length=Nx)), digits=15)
        ys = round.(collect(range(miny, stop=maxy, length=Ny)), digits=15)

        pos = zeros(Float64, 2, Nx, Ny)
        for i in 1:Nx, j in 1:Ny
            pos[:, i, j] .= (xs[i], ys[j])
        end

        v = zeros(Float64, (2,Nx,Ny))
        v_new = zeros(Float64, (2,Nx,Ny))
        momentum = zeros(Float64, (2,Nx,Ny))
        momentum_new = zeros(Float64, (2,Nx,Ny))
        f_ext = zeros(Float64, (2,Nx,Ny))
        f_int = zeros(Float64, (2,Nx,Ny))
        
        m = zeros(Float64, Nx, Ny)

        # empty_tuple_vector = Vector{Tuple{Int64, Int64}}()
    

        new(pos, 
        v, v_new, 
        momentum, momentum_new,
        f_ext, f_int,
        m, dx, dy, minx, miny)
    end
end


"""
MPMSimulation Type. Contains all necessary information about an MPM simulation.
"""
mutable struct MPMSimulation{MPGroupTouple <: Tuple}
    mp_groups::MPGroupTouple        # Tuple of MaterialPointGroups
    grid::Grid
    dt::Float64
    t::Float64
    total_time::Float64
end


