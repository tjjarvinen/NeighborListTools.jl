using BenchmarkTools
using NeighborListTools
using AtomsBuilder
using Unitful

SUITE = BenchmarkGroup()

ar_100_000 = bulk(:Ar, cubic=true) * 30
ar_500_000 = bulk(:Ar, cubic=true) * 50
ar_1_000_000 = bulk(:Ar, cubic=true) * 64

cutoff = 9.0u"Å"

SUITE["Cell list"] = BenchmarkGroup()
SUITE["Pair list"] = BenchmarkGroup()

SUITE["Cell list"]["100_000"] = @benchmarkable NeighborListTools.form_cell_list(
    $ar_100_000,
    $cutoff
)

SUITE["Cell list"]["500_000"] = @benchmarkable NeighborListTools.form_cell_list(
    $ar_500_000,
    $cutoff
)

SUITE["Cell list"]["1_000_000"] = @benchmarkable NeighborListTools.form_cell_list(
    $ar_1_000_000,
    $cutoff
)


cl_100_000 = NeighborListTools.form_cell_list(ar_100_000, cutoff)
cl_500_000 = NeighborListTools.form_cell_list(ar_500_000, cutoff)
cl_1_000_000 = NeighborListTools.form_cell_list(ar_1_000_000, cutoff)

SUITE["Pair list"]["100_000"] = @benchmarkable NeighborListTools.get_pairs(
    $cl_100_000,
    2,2,2
)

SUITE["Pair list"]["500_000"] = @benchmarkable NeighborListTools.get_pairs(
    $cl_500_000,
    2,2,2
)

SUITE["Pair list"]["1_000_000"] = @benchmarkable NeighborListTools.get_pairs(
    $cl_1_000_000,
    2,2,2
)