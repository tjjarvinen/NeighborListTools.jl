using BenchmarkTools
using NeighborListTools
using AtomsBuilder
using AtomsBase
using Unitful

using CellListMap
using NeighbourLists

##

ar_10_000 = bulk(:Ar, cubic=true) * 14
ar_55_000 = bulk(:Ar, cubic=true) * 25
ar_100_000 = bulk(:Ar, cubic=true) * 30
ar_500_000 = bulk(:Ar, cubic=true) * 50
ar_1_000_000 = bulk(:Ar, cubic=true) * 64

cutoff = 9.0u"Å"


##

function test_sitelist(sys, cutoff; use_fp32=false)
    cl = NeighborListTools.form_cell_list(sys, cutoff; use_fp32=use_fp32)
    s = size(cl)
    tmp = [ (i,j,k) for i in 1:s[1], j in 1:s[2], k in 1:s[2]]
    Threads.@threads for x in tmp
        NeighborListTools.get_site_pairlist(cl, x...)
    end
end

function test_pairlist(sys, cutoff; use_fp32=false)
    cl = NeighborListTools.form_cell_list(sys, cutoff; use_fp32=use_fp32)
    s = size(cl)
    tmp = [ (i,j,k) for i in 1:s[1], j in 1:s[2], k in 1:s[2]]
    Threads.@threads for x in tmp
        NeighborListTools.get_pairlist(cl, x...)
    end
end

function test_CellListMap(sys, cutoff)
    pos = position(sys, :)
    box = bounding_box(sys)
    bb = [box[1][1], box[2][2], box[3][3]]
    CellListMap.neighborlist(pos, cutoff; unitcell=bb)
end



##


SUITE = BenchmarkGroup()
SUITE["Cell list"] = BenchmarkGroup()
SUITE["Pairs for one cell"] = BenchmarkGroup()
SUITE["Pair interaction list"] = BenchmarkGroup()
SUITE["Site interaction list"] = BenchmarkGroup()
SUITE["Site interaction list"]["FP64"] = BenchmarkGroup()
SUITE["Site interaction list"]["FP32"] = BenchmarkGroup()

SUITE["CellListMap.jl"] = BenchmarkGroup()
SUITE["NeighbourLists.jl"] = BenchmarkGroup()

SUITE["Cell list"]["10_000"] = @benchmarkable NeighborListTools.form_cell_list(
    $ar_10_000,
    $cutoff
)

SUITE["Cell list"]["55_000"] = @benchmarkable NeighborListTools.form_cell_list(
    $ar_55_000,
    $cutoff
)


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

cl_10_000 = NeighborListTools.form_cell_list(ar_10_000, cutoff; use_fp32=false)
cl_55_000 = NeighborListTools.form_cell_list(ar_55_000, cutoff; use_fp32=false)
cl_100_000 = NeighborListTools.form_cell_list(ar_100_000, cutoff; use_fp32=false)
cl_500_000 = NeighborListTools.form_cell_list(ar_500_000, cutoff; use_fp32=false)
cl_1_000_000 = NeighborListTools.form_cell_list(ar_1_000_000, cutoff; use_fp32=false)

SUITE["Pairs for one cell"]["center"] = @benchmarkable NeighborListTools.get_site_pairlist(
    $cl_100_000,
    2,2,2
)

SUITE["Pairs for one cell"]["plane surface"] = @benchmarkable NeighborListTools.get_site_pairlist(
    $cl_500_000,
    2,2,1
)

SUITE["Pairs for one cell"]["two plane intersection"] = @benchmarkable NeighborListTools.get_site_pairlist(
    $cl_1_000_000,
    2,1,1
)

SUITE["Pairs for one cell"]["corner"] = @benchmarkable NeighborListTools.get_site_pairlist(
    $cl_1_000_000,
    1,1,1
)

SUITE["Site interaction list"]["FP64"]["10_000"] = @benchmarkable test_sitelist($ar_10_000, $cutoff)
SUITE["Site interaction list"]["FP64"]["55_000"] = @benchmarkable test_sitelist($ar_55_000, $cutoff)
SUITE["Site interaction list"]["FP64"]["100_000"] = @benchmarkable test_sitelist($ar_100_000, $cutoff)
SUITE["Site interaction list"]["FP64"]["500_000"] = @benchmarkable test_sitelist($ar_500_000, $cutoff)
SUITE["Site interaction list"]["FP64"]["1_000_000"] = @benchmarkable test_sitelist($ar_1_000_000, $cutoff)

SUITE["Site interaction list"]["FP32"]["10_000"] = @benchmarkable test_sitelist($ar_10_000, $cutoff; use_fp32=$true)
SUITE["Site interaction list"]["FP32"]["55_000"] = @benchmarkable test_sitelist($ar_55_000, $cutoff; use_fp32=$true)
SUITE["Site interaction list"]["FP32"]["100_000"] = @benchmarkable test_sitelist($ar_100_000, $cutoff; use_fp32=$true)
SUITE["Site interaction list"]["FP32"]["500_000"] = @benchmarkable test_sitelist($ar_500_000, $cutoff; use_fp32=$true)
SUITE["Site interaction list"]["FP32"]["1_000_000"] = @benchmarkable test_sitelist($ar_1_000_000, $cutoff; use_fp32=$true)

SUITE["Pair interaction list"]["10_000"] = @benchmarkable test_pairlist($ar_10_000, $cutoff)
SUITE["Pair interaction list"]["55_000"] = @benchmarkable test_pairlist($ar_55_000, $cutoff)
SUITE["Pair interaction list"]["100_000"] = @benchmarkable test_pairlist($ar_100_000, $cutoff)
SUITE["Pair interaction list"]["500_000"] = @benchmarkable test_pairlist($ar_500_000, $cutoff)
SUITE["Pair interaction list"]["1_000_000"] = @benchmarkable test_pairlist($ar_1_000_000, $cutoff)

SUITE["NeighbourLists.jl"]["10_000"] = @benchmarkable NeighbourLists.PairList($ar_10_000, $cutoff)
SUITE["NeighbourLists.jl"]["55_000"] = @benchmarkable NeighbourLists.PairList($ar_55_000, $cutoff)
SUITE["NeighbourLists.jl"]["100_000"] = @benchmarkable NeighbourLists.PairList($ar_100_000, $cutoff)
SUITE["NeighbourLists.jl"]["500_000"] = @benchmarkable NeighbourLists.PairList($ar_500_000, $cutoff)
SUITE["NeighbourLists.jl"]["1_000_000"] = @benchmarkable NeighbourLists.PairList($ar_1_000_000, $cutoff)

SUITE["CellListMap.jl"]["10_000"] = @benchmarkable test_CellListMap($ar_10_000, $cutoff)
SUITE["CellListMap.jl"]["55_000"] = @benchmarkable test_CellListMap($ar_55_000, $cutoff)
SUITE["CellListMap.jl"]["100_000"] = @benchmarkable test_CellListMap($ar_100_000, $cutoff)
SUITE["CellListMap.jl"]["500_000"] = @benchmarkable test_CellListMap($ar_500_000, $cutoff)
SUITE["CellListMap.jl"]["1_000_000"] = @benchmarkable test_CellListMap($ar_1_000_000, $cutoff)