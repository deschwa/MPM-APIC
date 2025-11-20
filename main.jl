include("src/mpm.jl")
using .MPM
using Base.Threads

include("scripts/read_csv.jl")
include("scripts/write_csv.jl")

input_csv = "input/input.csv"
input_yaml = "input/input.yaml"


particle_path(t) = "output/particles/dump_p.$t.xyz"
grid_path(t) = "output/grid/dump_g.$t.xyz"


sim, grid_dict, mat_dict = create_sim_from_csv(input_csv, input_yaml)


timestep!(sim, 0.0)

println("Starting simulation with $(Threads.nthreads()) Thread(s)...")
calculated_steps = 0
start_time = time()
while sim.t < sim.total_time
    timestep!(sim, 0.1)
    # print("Calculating time: ", round(sim.t, digits=4), " or $(round(sim.t / sim.total_time, digits=4)*100)%\r")

    if calculated_steps % 50 == 0
        write_particle_xyz(sim, particle_path(calculated_steps))
        write_grid_xyz(sim, grid_path(calculated_steps))
    end

    global calculated_steps += 1
end
println("Simulation completed in $(round(time()-start_time, digits=4))s.                             ")
# zip_folder("output/particles", "output/particles.zip")
# zip_folder("output/grid", "output/grid.zip")

# write_particle_csv(sim, mat_dict, output_particles)
# write_grid_csv(sim, output_grid)