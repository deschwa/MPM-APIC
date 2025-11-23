include("src/mpm.jl")
using .MPM
using Base.Threads

using BenchmarkTools
using Profile

include("scripts/read_csv.jl")
include("scripts/write_csv.jl")

input_csv = "input/input.csv"
input_yaml = "input/input.yaml"


particle_path(t) = "output/particles/dump_p.$t.xyz"
grid_path(t) = "output/grid/dump_g.$t.xyz"


sim, grid_dict, mat_dict = create_sim_from_csv(input_csv, input_yaml)


timestep!(sim, 0.99)

function alt_shape_function(r_rel::SVector{2,Float64}, dx::Float64, dy::Float64)
    N_x = abs(r_rel[1]) <= dx ? 0.5 * (1.0 + cos(pi * r_rel[1] / dx)) : 0.0
    N_y = abs(r_rel[2]) <= dy ? 0.5 * (1.0 + cos(pi * r_rel[2] / dy)) : 0.0
    N_I = N_x * N_y

    dNxdx = abs(r_rel[1]) <= dx ? 0.5 * pi * sin(pi * r_rel[1] / dx) / dx : 0.0
    dNydy = abs(r_rel[2]) <= dy ? 0.5 * pi * sin(pi * r_rel[2] / dy) / dy : 0.0
    ∇N_I = @SVector [dNxdx * N_y,
            dNydy * N_x]
    return (N_I, ∇N_I)
       
end

alt_shape_function(SVector(0.05, 0.05), 0.1, 0.1)



println("Benchmarking original shape function:")
display(@benchmark shape_function(r, 0.1, 0.1) setup=(r = @SVector rand(2)) samples=100000 evals=1)


println("Benchmarking alternative shape function:")
display(@benchmark alt_shape_function(r, 0.1, 0.1) setup=(r = @SVector rand(2)) samples=100000 evals=1)