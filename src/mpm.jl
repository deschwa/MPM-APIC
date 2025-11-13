module MPM

include("types.jl")
include("shapefunctions.jl")
include("helpers.jl")
include("materials.jl")
include("boundaries.jl")
include("timestep.jl")


export MPMSimulation, MaterialPointGroup, Grid, timestep!


end # Nodule MPM