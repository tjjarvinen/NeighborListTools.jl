using BenchmarkTools
using NeighborListTools
using AtomsBuilder
using Unitful

##

ar_100_000 = bulk(:Ar, cubic=true) * 30
ar_500_000 = bulk(:Ar, cubic=true) * 50
ar_1_000_000 = bulk(:Ar, cubic=true) * 64

cutoff = 9.0u"Å"


##

function test_sitelist(cl)
    s = size(cl)
    tmp = [ (i,j,k) for i in 1:s[1], j in 1:s[2], k in 1:s[2]]
    Threads.@threads for x in tmp
        NeighborListTools.get_pairs(cl, x...)
    end
end

function test_pairlist(cl)
    s = size(cl)
    tmp = [ (i,j,k) for i in 1:s[1], j in 1:s[2], k in 1:s[2]]
    Threads.@threads for x in tmp
        NeighborListTools.get_pair_list(cl, x...)
    end
end



##


SUITE = BenchmarkGroup()
SUITE["Cell list"] = BenchmarkGroup()
SUITE["Pairs for one cell"] = BenchmarkGroup()
SUITE["Pair interaction list"] = BenchmarkGroup()
SUITE["Site interaction list"] = BenchmarkGroup()

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

SUITE["Pairs for one cell"]["center"] = @benchmarkable NeighborListTools.get_pairs(
    $cl_100_000,
    2,2,2
)

SUITE["Pairs for one cell"]["plane surface"] = @benchmarkable NeighborListTools.get_pairs(
    $cl_500_000,
    2,2,1
)

SUITE["Pairs for one cell"]["two plane intersection"] = @benchmarkable NeighborListTools.get_pairs(
    $cl_1_000_000,
    2,1,1
)

SUITE["Pairs for one cell"]["corner"] = @benchmarkable NeighborListTools.get_pairs(
    $cl_1_000_000,
    1,1,1
)


SUITE["Site interaction list"]["100_000"] = @benchmarkable test_sitelist($cl_100_000)
SUITE["Site interaction list"]["500_000"] = @benchmarkable test_sitelist($cl_500_000)
SUITE["Site interaction list"]["1_000_000"] = @benchmarkable test_sitelist($cl_1_000_000)

SUITE["Pair interaction list"]["100_000"] = @benchmarkable test_pairlist($cl_100_000)
SUITE["Pair interaction list"]["500_000"] = @benchmarkable test_pairlist($cl_500_000)
SUITE["Pair interaction list"]["1_000_000"] = @benchmarkable test_pairlist($cl_1_000_000)