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

