# include("../src/mpm.jl")
# using .MPM


using CSV
using DataFrames
using StaticArrays
using YAML


function write_particle_csv(sim::MPMSimulation, material_dict::Dict{String, Any}, output_path::String)
    rows = []
    for mp_group in sim.mp_groups
        for i in 1:length(mp_group.mass)
            row = (
                type = mp_group.type,
                x = mp_group.pos[i][1],
                y = mp_group.pos[i][2],
                vx = mp_group.vel[i][1],
                vy = mp_group.vel[i][2],
                mass = mp_group.mass[i],
                volume = mp_group.volume[i],
            )
            push!(rows, row)
        end
    end
    df = DataFrame(rows)
    CSV.write(output_path, df)

end



function write_grid_csv(sim::MPMSimulation, output_path::String)
    rows = []
    grid = sim.grid
    Nx, Ny = size(grid.pos)
    for i in 1:Nx
        for j in 1:Ny
            row = (
                x = grid.pos[i,j][1],
                y = grid.pos[i,j][2],
                vx = grid.v[i,j][1],
                vy = grid.v[i,j][2],
                mass = grid.mass[i,j],
            )
            push!(rows, row)
        end
    end
    df = DataFrame(rows)
    CSV.write(output_path, df)
end