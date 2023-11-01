using AtomsIO
using AtomsBase
using NeighbourLists
using NeighborListTools
using Unitful

##
fname_xyz = joinpath(pkgdir(NeighborListTools), "data", "TiAl-big.xyz")

# Example data look https://github.com/JuliaMolSim/AtomsBase.jl for details
# Here we convert is to FastSystem that has better data structure for our case
data = load_system(fname_xyz) |> FastSystem

# Main features
#
# atom positions
r = position(data)

# Cell vectors
cell = bounding_box(data)

# And individual atoms with index
a1 = data[1]

## Reference neighborlist

# this function uses existing package to calculate reference neighborlist
function neighborlist(ab, cutoff; length_unit=u"Å", kwargs...)
    cell = ustrip.(length_unit, hcat( bounding_box(ab)... )' )
    pbc = map( boundary_conditions(ab) ) do x
        x == Periodic()
    end
    r = map( 1:length(ab)) do i
        # Need to have SVector here for PairList to work
        # Here we use FastSystem to get it.
        ustrip.(length_unit, position(ab,i))
    end
    nlist = PairList(r, ustrip(length_unit, cutoff), cell, pbc; int_type=Int)
    return nlist
end


# Here is an example neighborlist - look ?PairList for more info
nlist = neighborlist(data, 5.5u"Å")

# Distance vectors
nlist.X

# Atom indices in original data
nlist.i, nlist,j

# This is one possible output, which we should implement, but we should also make other types of output.