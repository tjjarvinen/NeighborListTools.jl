

function get_pairs(cl::CellList, indx1::CartesianIndex, indx2::CartesianIndex)
    c1 = cl[indx1]
    c2 = cl[indx2]
    return _get_sites(cl, c1, c2)
end

function get_site_pairlist(cl::CellList, i::Int, j::Int, k::Int)
    c1 = cl[i,j,k]
    c2 = cl[i-1:i+1, j-1:j+1, k-1:k+1]
    return _get_sites(cl, c1, c2)
end

function get_site_pairlist(cl::CellList, carindex::CartesianIndex)
    return get_site_pairlist(cl, carindex[1], carindex[2], carindex[3])
end

function _get_sites(cl::CellList, c1, c2)
    fi1, fi2 = _distances_sites(cl.cutoff, c1.positions, c2.positions)

    return (
        indx1 = fi1,
        indx2 = fi2,
        origin_indx1 = view(c1.indices, fi1),
        origin_indx2 = view(c2.indices, fi2),
        species1 = view(c1.species,fi1),
        species2 = view(c2.species, fi2),
        r = view(c2.positions,fi2) - view(c1.positions,fi1),
    )
end

function _get_pairs(cl::CellList, c1, c2)
    fi1, fi2 = _distances_pairs(cl.cutoff, c1.positions, c2.positions)

    return (
        indx1 = fi1,
        indx2 = fi2,

        #origin_indx1 = c1.indices[fi1],
        #origin_indx2 = c2.indices[fi2],
        #species1 = c1.species[fi1],
        #species2 = c2.species[fi2],
        #r = view(c2.positions,fi2) - view(c1.positions,fi1)

        # This is a much faster, but is it optimal for the later use?
        origin_indx1 = view(c1.indices, fi1),
        origin_indx2 = view(c2.indices, fi2),
        species1 = view(c1.species,fi1),
        species2 = view(c2.species, fi2),
        r = view(c2.positions,fi2) - view(c1.positions,fi1),
    )
end

function get_pairlist(cl::CellList, i::Int, j::Int, k::Int)
    c1 = cl[i,j,k]
    c2 = cl[i:i+1, j:j+1, k:k+1]
    return _get_pairs(cl, c1, c2)
end

function get_pairlist(cl::CellList, carindex::CartesianIndex)
    return get_pairlist(cl, carindex[1], carindex[2], carindex[3])
end


## Kernels for pair index generation

# generate indeces that are then checked
@kernel function make_idx_kernel(out, @Const(R))
    i, j = @index(Global, NTuple)
    @inbounds out[i,j] = (i,j)
end

function make_idx(idx)
    backend = get_backend(idx)
    out = allocate(backend, Tuple{Int32, Int32}, size(idx)...)
    kernel = make_idx_kernel(backend)
    kernel(out, idx; ndrange = size(out))
    return out[idx]
end

# test for conditions
function site_test(R, cutoff)
    out = @. zero(cutoff) < R < cutoff
    return out
end


@kernel function pair_test_kernel(out, @Const(R), rcut)
    i, j = @index(Global, NTuple)
    @inbounds out[i,j] = R[i,j] < rcut && i < j
end

function pair_test(R, cutoff)
    backend = get_backend(R)
    out = similar(R, Bool)

    kernel = pair_test_kernel(backend)
    kernel(out, R, cutoff; ndrange = size(R))
    return out
end

function _dis(r1, r2)
    Δr = r2 - r1
    return sqrt( sum(x->x^2, Δr) )
end


@kernel function _distances_kernel(out, @Const(R1), @Const(R2))
    J = @index(Global)
    for i in 1:length(R1)
        @inbounds out[i, J] = _dis(R1[i], R2[J])
    end
end

# General CPU + GPU
function _distances_pairs(cutoff, r1, r2)
    backend = get_backend(r1)
    @assert backend == get_backend(r2)

    d = similar(r1, (eltype∘eltype)(r1), length(r1), length(r2))

    kernel = _distances_kernel(backend)
    kernel(d, r1, r2; ndrange=length(r2))
    #synchronize(backend)
    test_results = pair_test(d, cutoff)
    idx = make_idx(test_results)
    i = map(x->x[1], idx)
    j = map(x->x[2], idx)
    return i, j
end

function _distances_sites(cutoff, r1, r2)
    backend = get_backend(r1)
    @assert backend == get_backend(r2)

    d = similar(r1, (eltype∘eltype)(r1), length(r1), length(r2))

    kernel = _distances_kernel(backend)
    kernel(d, r1, r2; ndrange=length(r2))
    #synchronize(backend)
    test_results = site_test(d, cutoff)
    idx = make_idx(test_results)
    i = map(x->x[1], idx)
    j = map(x->x[2], idx)
    return i, j
end

# For CPU
function _distances_pairs(
    cutoff, 
    r1::Union{SubArray{<:Any, 1, <:Array, <:Any}, Vector}, 
    r2::Union{SubArray{<:Any, 1, <:Array, <:Any}, Vector}
)
    r = Distances.pairwise(Distances.Euclidean(), r1, r2)
    findx = findall(x-> x <(cutoff), r)
    fi1 = [ x[1] for x in findx if x[1] < x[2] ]
    fi2 = [ x[2] for x in findx if x[1] < x[2] ]
    return fi1, fi2
end

function _distances_sites(
    cutoff,
    r1::Union{SubArray{<:Any, 1, <:Array, <:Any}, Vector}, 
    r2::Union{SubArray{<:Any, 1, <:Array, <:Any}, Vector}
)
    r = Distances.pairwise(Distances.Euclidean(), r1, r2)
    findx = findall(x-> zero(cutoff) < x <(cutoff), r)
    fi1 = [ x[1] for x in findx ]
    fi2 = [ x[2] for x in findx ]
    return fi1, fi2
end

## end kernels for pair index generation