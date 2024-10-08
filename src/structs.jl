

mutable struct CellList{TC, TP, TS, TI, TCO, TCE}
    cell::TCE
    cell_indices::TC
    cutoff::TCO
    indices::TI
    positions::TP
    species::TS
end

# Change array type for CellList
function CellList(TA, cl::CellList)
    ta_positions = TA(cl.positions)
    ta_species = TA(cl.species)
    ta_indices = TA(cl.indices)
    return CellList(
        cl.cell,
        cl.cell_indices,
        cl.cutoff,
        ta_indices,
        ta_positions,
        ta_species
    )
end

Base.size(cl::CellList) = size(cl.cell_indices)
Base.length(cl::CellList) = length(cl.cell_indices)

Base.axes(cl::CellList, i) = axes(cl.cell_indices, i)
Base.IndexStyle(::CellList) = IndexCartesian()
Base.keys(cl::CellList) = CartesianIndices(size(cl))



function Base.getindex(cl::CellList, carindex::CartesianIndex)
    return cl[carindex[1], carindex[2], carindex[3]]    
end

function Base.getindex(cl::CellList, i)
    indx = vcat(i) do j
        cl[j]
    end
    return (
        positions = view(cl.positions, indx),
        species = view(cl.species, indx),
        indices = view(cl.indices, indx)
    )
end

function Base.getindex(cl::CellList, i::Int)
    indx = cl.cell_indices[i]
    return (
        positions = view(cl.positions, indx),
        species = view(cl.species, indx),
        indices = view(cl.indices, indx)
    )
end

function Base.getindex(cl::CellList, i::Int, j::Int, k::Int)
    s = size(cl)
    if 0<i<=s[1] && 0<j<=s[2] && 0<k<=s[3]
        indx = cl.cell_indices[i,j,k]
        return (
            positions = view(cl.positions, indx),
            species = view(cl.species, indx),
            indices = view(cl.indices, indx),
            sift = _get_shift(cl, i,j,k)
        )
    else
        ii = i>0 ? rem(i-1, s[1]) + 1 : i%s[1] + s[1]
        jj = j>0 ? rem(j-1, s[2]) + 1 : j%s[2] + s[2]
        kk = k>0 ? rem(k-1, s[3]) + 1 : k%s[3] + s[3]
        indx = cl.cell_indices[ii,jj,kk]
        sift = _get_shift(cl, i,j,k)
        return (
            positions = view(cl.positions, indx),
            species = view(cl.species, indx),
            indices = view(cl.indices, indx),
            sift = sift
        )
    end
end

function _get_shift(cl, i::Int, j::Int, k::Int)
    s = size(cl)
    ii = i>0 ? div(i-1, s[1]) : div(i,s[1]) -1
    jj = j>0 ? div(j-1, s[2]) : div(j,s[2]) -1
    kk = k>0 ? div(k-1, s[3]) : div(k,s[3]) -1
    Δr = ii*cl.cell[1] + jj*cl.cell[2] + kk*cl.cell[3]
    return Δr
end


function Base.getindex(cl::CellList, i, j, k)
    s = size(cl)
    if 1<=minimum(i) && maximum(i)<=s[1] &&
            1<=minimum(j) && maximum(j)<=s[2] &&
            1<=minimum(k) && maximum(k)<=s[3]
        indx = reduce(vcat, cl.cell_indices[i,j,k])
        return (
            positions = view(cl.positions, indx),
            species = view(cl.species, indx),
            indices = view(cl.indices, indx)
        )
    else
        tmp = [ cl[a,b,c] for a in i, b in j, c in k ]
        pos = mapreduce(vcat, tmp) do x
            map(y->y - x.sift, x.positions)
        end
        spc = mapreduce(x->x.species, vcat, tmp)
        indx = mapreduce(x->x.indices, vcat, tmp)
        return (
            positions = pos,
            species = spc,
            indices = indx
        )
    end
end