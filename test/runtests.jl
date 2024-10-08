using AtomsBuilder
using NeighborListTools
using Unitful
using Test

##


@testset "basics" begin
    cutoff = 9.0u"Å"

    sys = bulk(:Ar, cubic=true) * 50
    q = NeighborListTools.form_cell_list(sys, cutoff)
    w = q[1,1,1]

    @test haskey(w, :positions)
    @test haskey(w, :indices)
    @test haskey(w, :species)
    @test haskey(w, :sift)

    pairs = NeighborListTools.get_site_pair_list(q, 1,1,1)

    @test haskey(pairs, :r)
    @test haskey(pairs, :species1)
    @test haskey(pairs, :species2)
    @test haskey(pairs, :origin_indx1)
    @test haskey(pairs, :origin_indx2)
    @test haskey(pairs, :indx1)
    @test haskey(pairs, :indx2)

    pairs = NeighborListTools.get_pair_list(q, 1,1,1)

    @test haskey(pairs, :r)
    @test haskey(pairs, :species1)
    @test haskey(pairs, :species2)
    @test haskey(pairs, :origin_indx1)
    @test haskey(pairs, :origin_indx2)
    @test haskey(pairs, :indx1)
    @test haskey(pairs, :indx2)
end