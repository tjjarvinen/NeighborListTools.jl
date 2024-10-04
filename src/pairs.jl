

function get_pairs(cl::CellList, i::Int, j::Int, k::Int)
    c1 = cl[i,j,k]
    c2 = cl[[i-1,i+1], [j-1,j+1], [k-1,k+1]]

    r = pairwise(Euclidean(), c1.positions, c2.positions)
    findx = findall(<(cl.cutoff), r)
    fi1 = [ x[1] for x in findx ]
    fi2 = [ x[2] for x in findx ]

    r2 = pairwise(Euclidean(), c1.positions, c1.positions)
    findx2 = findall(x -> zero(cl.cutoff) < x < cl.cutoff, r2)
    ff1 = [ x[1] for x in findx2 ]
    ff2 = [ x[1] for x in findx2 ]

    return (
        distances = vcat(r[findx], r2[findx2]),
        ind1 = vcat(c1.indices[fi1],  c1.indices[ff1]),
        ind2 = vcat(c2.indices[fi2], c1.indices[ff2]),
        r = vcat(
            c2.positions[fi2] - c1.positions[fi1],
            c1.positions[ff2] - c1.positions[ff1]
        )
    )
end


function get_pairs(cl::CellList, indx1, indx2)
    c1 = cl[indx1]
    c2 = cl[indx2]

    r = pairwise(Euclidean(), c1.positions, c2.positions)
    findx = findall(x-> zero(cl.cutoff)< x <(cl.cutoff), r)
    fi1 = [ x[1] for x in findx ]
    fi2 = [ x[2] for x in findx ]

    return (
        distances = r[findx],
        ind1 = c1.indices[fi1],
        ind2 = c2.indices[fi2],
        r = c2.positions[fi2] - c1.positions[fi1],
    )
end