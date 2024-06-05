using AtomsIO, AtomsBase, Unitful
using NeighborListTools
using LinearAlgebra, StaticArrays
using Test

# TODO: Turn into @testset style tests 

# Example data look https://github.com/JuliaMolSim/AtomsBase.jl for details
# Here we convert is to FastSystem that has better data structure for our case
fname_xyz = joinpath(pkgdir(NeighborListTools), "data", "TiAl-big.xyz")
data = load_system(fname_xyz) |> FastSystem

# user inputs
cutoff = 4.0u"Å" # example cutoff
cell_size_suggestion = [4.2u"Å", 4.3u"Å", 4.0u"Å"] # example cell size suggestions. must be bigger than cutoff in each dimension

# Test unsorted CellList object
c = CellList(data, cutoff, cell_size_suggestion)
# inspect that they were not sorted
c.celllist[1].position
c.celllist[1].elementtypes
c.celllist[1].originalindex
c.celllist[1].cellid

# original atoms are in what cell?
c.originalindstocell # We see atom #1 is now in cell (1,2,2)
getfield.(c.celllist,:cellid)
c.celllist[c.originalindstocell[1]...].originalindex # and we find that atom in (1,2,2) as expected

# Test sorting an unsorted CellList
c.celllist[1].sorttype
cellsort!.(c.celllist, :atomicnumber)
c.celllist[1].sorttype
c.celllist[1].position
c.celllist[1].elementtypes
c.celllist[1].originalindex
c.originalindstocell

# Test generating a sorted CellList directly
c2 = CellList(data, cutoff, cell_size_suggestion, sortstyle=:alphaelement)
c2.celllist[1].sorttype == c.celllist[1].sorttype
c2.celllist[1].position == c.celllist[1].position
c2.celllist[1].elementtypes == c.celllist[1].elementtypes
c2.celllist[1].originalindex == c.celllist[1].originalindex
c2.originalindstocell
natoms(c) == natoms(c2)


using BenchmarkTools, ProfileView
natoms(c)
@btime CellList($data, $cutoff, $cell_size_suggestion)
ProfileView.@profview for j = 1:10000
    CellList(data, cutoff, cell_size_suggestion)
end