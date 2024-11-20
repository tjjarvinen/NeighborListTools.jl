using AtomsBuilder
using NeighborListTools
using Unitful
using Test

##


@testset "basics" begin
    cutoff = 9.0u"Å"

    sys = bulk(:Ar, cubic=true) * 20
    q = CellList(sys, cutoff)
    w = q[1,1,1]

    @test haskey(w, :positions)
    @test haskey(w, :indices)
    @test haskey(w, :species)
    @test haskey(w, :sift)

    pairs = get_site_pairlist(q, 1,1,1)

    @test haskey(pairs, :r)
    @test haskey(pairs, :species1)
    @test haskey(pairs, :species2)
    @test haskey(pairs, :origin_indx1)
    @test haskey(pairs, :origin_indx2)
    @test haskey(pairs, :indx1)
    @test haskey(pairs, :indx2)

    pairs = get_pairlist(q, 1,1,1)

    @test haskey(pairs, :r)
    @test haskey(pairs, :species1)
    @test haskey(pairs, :species2)
    @test haskey(pairs, :origin_indx1)
    @test haskey(pairs, :origin_indx2)
    @test haskey(pairs, :indx1)
    @test haskey(pairs, :indx2)
end


@testset "Iterator interface" begin
    cutoff = 9.0u"Å"
    sys = bulk(:Ar, cubic=true) * 20
    cl = CellList(sys, cutoff)

    @testset "Pair iterator"  begin
        pair_iter = give_pair_iterators(cl, 4)
        
        @test length(pair_iter) == 4

        new, i = iterate(pair_iter[1])

        @test isa(new, NamedTuple)
        @test hasproperty(new, :r)
        @test hasproperty(new, :species1)
        @test hasproperty(new, :species2)
        @test hasproperty(new, :origin_indx1)
        @test hasproperty(new, :origin_indx2)

        @test isa(new.r, AbstractVector)    
        @test isa(new.species1, AbstractVector)
        @test isa(new.species2, AbstractVector)
        @test isa(new.origin_indx1, AbstractVector)
        @test isa(new.origin_indx2, AbstractVector)

        @test length(new.r) == length(new.species1)
        @test length(new.r) == length(new.species2)
        @test length(new.r) == length(new.origin_indx1)
        @test length(new.r) == length(new.origin_indx2)
    end

    @testset "Site iterator" begin
        site_iter = give_site_iterators(cl, 4)

        @test length(site_iter) == 4

        new, i = iterate(site_iter[1])

        @test isa(new, NamedTuple)
        @test hasproperty(new, :r)
        @test hasproperty(new, :species1)
        @test hasproperty(new, :species2)
        @test hasproperty(new, :origin_indx1)
        @test hasproperty(new, :origin_indx2)

        @test isa(new.r, AbstractVector)    
        @test isa(new.species1, AbstractVector)
        @test isa(new.species2, AbstractVector)
        @test isa(new.origin_indx1, AbstractVector)
        @test isa(new.origin_indx2, AbstractVector)

        @test length(new.r) == length(new.species1)
        @test length(new.r) == length(new.species2)
        @test length(new.r) == length(new.origin_indx1)
        @test length(new.r) == length(new.origin_indx2)
    end
    
end