using LinearAlgebra
include("src/mpm.jl")
using .MPM
using Base.Threads
using Plots

include("scripts/write_csv.jl")

particle_path(t) = "output/particles/dump_p.$(lpad(string(t), 5, '0')).xyz"



function rotating_square(midpoint::AbstractArray, side_length::Float64, N_1d::Int64, angular_speed::Float64)
    # Initialize
    positions = zeros(Float64, 2, N_1d*N_1d)
    velocities = zeros(Float64, 2, N_1d*N_1d)

    minx = midpoint[1] - 0.5 * side_length
    maxx = midpoint[1] + 0.5 * side_length
    miny = midpoint[2] - 0.5 * side_length
    maxy = midpoint[2] + 0.5 * side_length

    dx = side_length/N_1d
    volume = dx^2
    volumes = [volume for _ in 1:N_1d^2]

    xs = collect(range(minx, maxx, length=N_1d))
    ys = collect(range(miny, maxy, length=N_1d))

    # @assert length(xs) == size(positions[2]) "Shapes do not match!"
    @assert length(xs) == length(ys) "xs != ys"

    for i in 1:N_1d, j in 1:N_1d
        positions[:, i + (j-1)*N_1d] .= (xs[i], ys[j])
    end

    # Velocities
    for i in 1:length(positions[1,:])
        rel_vec = positions[:, i] - midpoint
        rel_vec_3d = [rel_vec..., 0.0]
        r = norm(rel_vec)
        v = r*angular_speed
        e_phi = cross([0.0, 0.0, 1.0], rel_vec_3d)[1:2]
        e_phi ./= norm(e_phi)
        v_vec = v * e_phi
        velocities[:,i] .= v_vec

        
    end

    return positions, velocities, volumes

end

function total_angular_momentum(mp_group::MaterialPointGroup)
    ang_mom = 0.0
    for i in 1:mp_group.N
        ang_mom += norm(cross([mp_group.pos[:,i]..., 0.0], [mp_group.vel[:,i]..., 0.0]))
    end
    return ang_mom
end


pos, vel, vols = rotating_square([0.0,0.0], 1.0, 10, 1.0)

mass = [1.0 for _ in 1:10^2]

LinEla = LinearElastic(1e3, 0.3, 1000)


mp = MaterialPointGroup(pos, vel, mass, vols, LinEla, "1")
grid = Grid(10, 10, -1.0, 1.0, -1.0, 1.0)

sim_template = MPMSimulation((mp,), grid, 1e-3, 0.0, 10.0)

alphas = collect(range(0.95, 1.0, length=11))

angular_momenta = Dict{Real, AbstractVector}()
times = Dict{Real, AbstractVector}()

plt = plot(title="Drehimpuls Erhaltung APIC", xlabel="t", ylabel="L")
for alpha in alphas
    sim = deepcopy(sim_template)
    calculated_steps = 0
    angular_momenta[alpha] = []
    times[alpha] = []

    while sim.t < sim.total_time
        timestep!(sim, alpha)
        if calculated_steps % 10 == 0
            print("Calculating time: ", round(sim.t, digits=1), " or $(round(sim.t / sim.total_time *100, digits=1))%\r")
        end

        # if calculated_steps % 50 == 0
        #     # println(particle_path(calculated_steps))
            
        #     write_particle_xyz(sim, particle_path(calculated_steps))
        # end
        L_current = total_angular_momentum(sim.mp_groups[1])
        push!(angular_momenta[alpha], L_current)
        push!(times[alpha], sim.t)
        calculated_steps += 1
    end
    # println("$(angular_momenta[alpha][1]), $(angular_momenta[alpha][end])")

    plot!(plt, times[alpha], angular_momenta[alpha], label="alpha=$alpha", legend=true)

    println("Plotting completed for alpha=$alpha")

end

# xlabel!("t")
# ylabel!("angular momentum")
# ylims!(0.0, 25.0)

savefig("angular_momenta.png")