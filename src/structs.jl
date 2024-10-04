

mutable struct CellList{TC, TP, TS, TI, TCO, TCE}
    cell::TCE
    cell_indices::TC
    cutoff::TCO
    indices::TI
    postions::TP
    species::TS
end

Base.size(cl::CellList) = size(cl.cell_indices)

Base.length(cl::CellList) = length(cl.cell_indices)


AtomsBase.cell(cl::CellList) = cl.cell

Base.IndexStyle(::CellList) = IndexCartesian()

Base.keys(cl::CellList) = CartesianIndices(size(cl))


function Base.getindex(cl::CellList, carindex::CartesianIndex)
    return cl[carindex[1], carindex[2], carindex[3]]    
end

function Base.getindex(cl::CellList, i)
    indx = vcat(cl.cell_indices[i]...)
    return (
        positions = view(cl.postions, indx),
        species = view(cl.species, indx),
        indices = view(cl.indices, indx)
    )
end

function Base.getindex(cl::CellList, i::Int)
    indx = cl.cell_indices[i]
    return (
        positions = view(cl.postions, indx),
        species = view(cl.species, indx),
        indices = view(cl.indices, indx)
    )
end

function Base.getindex(cl::CellList, i::Int, j::Int, k::Int)
    indx = cl.cell_indices[i,j,k]
    return (
        positions = view(cl.postions, indx),
        species = view(cl.species, indx),
        indices = view(cl.indices, indx)
    )
end

function Base.getindex(cl::CellList, i, j, k)
    indx = vcat(cl.cell_indices[i,j,k]...)
    return (
        positions = view(cl.postions, indx),
        species = view(cl.species, indx),
        indices = view(cl.indices, indx)
    )
end