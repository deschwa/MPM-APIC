# include("src/mpm.jl")
# using .MPM
# input_csv = "input/input.csv"
# input_yaml = "input/input.yaml"

using CSV
using DataFrames
using StaticArrays
using YAML


function get_sim_data(input_yaml)
    file = YAML.load_file(input_yaml)

    # Material dict
    mat_dict = Dict{String, Any}()
    for material in file["materials"]
        if material["type"] == "LinearElastic"
            mat_dict[material["name"]] = LinearElastic(
                material["E"], material["nu"], material["rho"]
            )
        elseif material["type"] == "NeoHookean"
            mat_dict[material["name"]] = NeoHookean(
                material["E"], material["nu"], material["rho"]
            )
        end
    end

    # Grid Data
    grid_dict = file["grid"][1]

    # time Data
    time_dict = file["time"][1]

    return mat_dict, grid_dict, time_dict
    
    
end

function create_particle_distribution_from_csv(input_file, material_dict, Nx, Ny)
    file = CSV.read(input_file, DataFrame)

    grouped_df = groupby(file, :type)

    mp_groups = []
    
    for (key, group) in pairs(grouped_df)
        type = key.type
        N = length(group.x)

        pos = zeros(Float64, 2, N)
        vel = zeros(Float64, 2, N)

        for i in 1:N
            pos[:,i] .= [Float64(group.x[i]), Float64(group.y[i])]
            vel[:,i] .= [Float64(group.vx[i]), Float64(group.vy[i])]
        end

        mass = [Float64(group.mass[i]) for i in 1:length(group.mass)]
        volume = [Float64(group.volume[i]) for i in 1:length(group.volume)]
        material = material_dict[string(type)]
        mp_group = MaterialPointGroup(pos, vel, mass, volume, material, string(type), Nx, Ny)
        push!(mp_groups, mp_group)
    end

    mp_groups = tuple(mp_groups...)
    return mp_groups
end


function create_grid_from_yaml(grid_dict)
    Nx = Int64(grid_dict["N_x"])
    Ny = Int64(grid_dict["N_y"])
    minx = Float64(grid_dict["min_x"])
    maxx = Float64(grid_dict["max_x"])
    miny = Float64(grid_dict["min_y"])
    maxy = Float64(grid_dict["max_y"])

    return Grid(Nx, Ny, minx, maxx, miny, maxy)
    
end


function create_sim_from_csv(input_csv, input_yaml)
    material_dict, grid_dict, time_dict = get_sim_data(input_yaml)
    grid = create_grid_from_yaml(grid_dict)
    Nx, Ny = grid_dict["N_x"], grid_dict["N_y"]
    mp_groups = create_particle_distribution_from_csv(input_csv, material_dict, Nx, Ny)
    dt = time_dict["dt"]
    t_end = time_dict["total_time"]

    sim = MPMSimulation(mp_groups, grid, dt, 0.0, t_end)
    return sim, grid_dict, material_dict
end


# create_sim_from_csv(input_csv, input_yaml)