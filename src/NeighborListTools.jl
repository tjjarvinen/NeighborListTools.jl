module NeighborListTools

using ArgCheck
using AtomsBase
using Distances
using KernelAbstractions
using LinearAlgebra
using StaticArrays
using Unitful

export Cell, CellList
export natoms, cellgrid, cellsort!, massivelyaddparticles!

include("structs.jl")
include("atoms_base_interface.jl")
include("pairs.jl")

# Old
#include("cell_list.jl")


end
