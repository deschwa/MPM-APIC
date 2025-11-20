# include("../src/mpm.jl")
# using .MPM


using CSV
using DataFrames
using StaticArrays
using YAML
using ZipFile
using Base.Filesystem

function write_particle_csv(sim::MPMSimulation, output_path::String)
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



function write_particle_xyz(sim::MPMSimulation, output_path::String)
    open(output_path, "w") do io
        total_particles = sum(length(mp_group.mass) for mp_group in sim.mp_groups)

        # Number of particles
        println(io, total_particles)

        # Comment line
        println(io, "type x y z vx vy vz σ_xx σ_yy σ_zz σ_xy σ_yz σ_zx volume")

        for mp_group in sim.mp_groups
            for i in 1:length(mp_group.mass)
                type = mp_group.type

                x = mp_group.pos[1, i]
                y = mp_group.pos[2, i]
                z = 0.0

                vx = mp_group.vel[1, i]
                vy = mp_group.vel[2, i]
                vz = 0.0

                σ_xx = mp_group.σ[1,1,i]
                σ_yy = mp_group.σ[2,2,i]
                σ_zz = 0.0
                σ_xy = mp_group.σ[1,2,i]
                σ_yz = 0.0
                σ_zx = 0.0

                volume = mp_group.volume[i]



                println(io, "$type $x $y $z $vx $vy $vz $σ_xx $σ_yy $σ_zz $σ_xy $σ_yz $σ_zx $volume")
            end
        end
    end
end


function write_grid_xyz(sim::MPMSimulation, output_path::String)
    grid = sim.grid
    Nx, Ny = size(grid.pos)
    open(output_path, "w") do io
        # Number of grid nodes
        println(io, Nx * Ny)

        # Comment line
        println(io, "x y z vx vy vz mass f_int_x f_int_y f_int_z f_ext_x f_ext_y f_ext_z")

        for i in 1:Nx
            for j in 1:Ny
                x = grid.pos[1, i,j]
                y = grid.pos[2, i,j]
                z = 0.0

                vx = grid.v[1, i,j]
                vy = grid.v[2, i,j]
                vz = 0.0

                mass = grid.mass[i,j]

                f_int_x = grid.f_int[1, i,j]
                f_int_y = grid.f_int[2, i,j]
                f_int_z = 0.0

                f_ext_x = grid.f_ext[1, i,j]
                f_ext_y = grid.f_ext[2, i,j]
                f_ext_z = 0.0
                
                println(io, "$x $y $z $vx $vy $vz $mass $f_int_x $f_int_y $f_int_z $f_ext_x $f_ext_y $f_ext_z")
            end
        end
    end
    
end



function zip_folder(source_dir::String, zip_file_name::String)
    @assert isdir(source_dir) "Das Quellverzeichnis existiert nicht: $source_dir"
    
    all_entries = readdir(source_dir)
    
    files_to_zip = filter(all_entries) do entry
        isfile(joinpath(source_dir, entry))
    end

    @assert !isempty(files_to_zip) "Keine Dateien zum Zippen im Verzeichnis: $source_dir"
 

    try
        ZipFile.writer(zip_file_name) do w
            for file in files_to_zip
                file_path = joinpath(source_dir, file)
                
                # Eine neue Datei im ZIP-Archiv erstellen
                # 'file' ist nur der Dateiname (z.B. "bild.jpg")
                f_in_zip = ZipFile.newfile(w, file)
                
                # Den Inhalt der Quelldatei in die ZIP-Datei schreiben
                write(f_in_zip, read(file_path))
                
                println("  - $file")
            end
        end
    catch e
        println("Fehler beim Erstellen der ZIP-Datei: $e")
        if isfile(zip_file_name)
            rm(zip_file_name)
            println("Unvollständige ZIP-Datei wurde entfernt.")
        end
    end
end