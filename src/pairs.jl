

function get_pairs(cl::CellList, indx1::CartesianIndex, indx2::CartesianIndex)
    c1 = cl[indx1]
    c2 = cl[indx2]
    return _get_pairs(cl, c1, c2)
end

function get_pairs(cl::CellList, i::Int, j::Int, k::Int)
    c1 = cl[i,j,k]
    c2 = cl[i-1:i+1, j-1:j+1, k-1:k+1]
    return _get_pairs(cl, c1, c2)
end


function _get_pairs(cl, c1, c2)
    r = pairwise(Euclidean(), c1.positions, c2.positions)
    findx = findall(x-> zero(cl.cutoff)< x <(cl.cutoff), r)
    fi1 = [ x[1] for x in findx ]
    fi2 = [ x[2] for x in findx ]

    return (
        distances = r[findx],
        origin_indx1 = c1.indices[fi1],
        origin_indx2 = c2.indices[fi2],
        indx1 = fi1,
        indx2 = fi2,
        r = c2.positions[fi2] - c1.positions[fi1],
    )
end

function get_pair_list(cl::CellList, i::Int, j::Int, k::Int)
    c1 = cl[i,j,k]
    c2 = cl[i:i+1, j:j+1, k:k+1]
    return _get_pairs(cl, c1, c2)
end