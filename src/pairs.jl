

function get_pairs(cl::CellList, indx1::CartesianIndex, indx2::CartesianIndex)
    c1 = cl[indx1]
    c2 = cl[indx2]
    return get_pairs(cl, c1, c2)
end

function get_site_pairlist(cl::CellList, i::Int, j::Int, k::Int)
    c1 = cl[i,j,k]
    c2 = cl[i-1:i+1, j-1:j+1, k-1:k+1]
    return get_pairs(cl, c1, c2)
end

function get_site_pairlist(cl::CellList, carindex::CartesianIndex)
    return get_site_pairlist(cl, carindex[1], carindex[2], carindex[3])
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
function _distances(r1, r2)
    backend = get_backend(r1)
    @assert backend == get_backend(r2)

    r = similar(r1, (eltype∘eltype)(r1), length(r1), length(r2))

    kernel = _distances_kernel(backend)
    kernel(r, r1, r2; ndrange=length(r2))
    return r
end

# For CPU
function _distances(r1::SubArray{<:Any, 1, <:Array, <:Any}, r2::SubArray{<:Any, 1, <:Array, <:Any})
    return Distances.pairwise(Distances.Euclidean(), r1, r2)
end


function get_pairs(cl::CellList, c1, c2)
    r = _distances(c1.positions, c2.positions)

    # This does not work with GPU
    findx = findall(x-> zero(cl.cutoff)< x <(cl.cutoff), r)
    
    fi1 = [ x[1] for x in findx ]
    fi2 = [ x[2] for x in findx ]

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
    return get_pairs(cl, c1, c2)
end

function get_pairlist(cl::CellList, carindex::CartesianIndex)
    return get_pairlist(cl, carindex[1], carindex[2], carindex[3])
end