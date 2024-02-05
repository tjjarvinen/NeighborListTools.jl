mutable struct Cell{TA, TI, TT}
    "positions of atoms"
    position::TA
    "indices in original input"
    originalindex::TI
    "how are the elements sorted? Currently only supports :alphaelement and :unsorted"
    sorttype::Symbol
    "vector of element types, stored as e.g. [Ti, Al, ...]"
    elementtypes::Vector{String}
    "assigned position of cell in a CellList if it was generated as part of one"
    cellid::TT
    # TODO: What other fields are needed?
end

# A method to create an empty cell with specific cellid
function Cell{T}(cellid::TT, sortstyle = :unsorted) where {T, TT}
    d = length(cellid)
    pos = Vector{SVector{d,T}}()
    originalindex = Vector{Int64}()
    elementranges = Vector{UnitRange{Int64}}()
    elementtypes = Vector{String}()
    return Cell{Vector{SVector{d,T}}, Vector{Int64}, TT}(pos, originalindex, sortstyle, elementtypes, cellid)
end

Base.show(io::IO, c::Cell) = print(io, "Cell with ", length(c), " atoms")
Base.show(io::IO, ::MIME"text/plain", c::Cell) = print(io, "Cell with ", length(c), " atoms")
Base.length(c::Cell) = length(c.position)

mutable struct CellList{T, TS, TB, TC, ND} # ND is number of dimensions
    "An array of Cell type objects"
    celllist::T
    "Axes and dimensions of the Cell objects in the CellList."
    cellshape::TS
    "Store the bounding box in AtomsBase compatible style."
    box::TB
    "Cutoff information used when constructing a neighbor list."
    cutoff::TC
    "Describe boundary conditions on each of the d sides of the bounding box."
    boundaryconditions::Symbol
    # TODO: What other fields are needed? Should the celllist also know the sorttype of all cells?
end

# TODO: There should be a way to get sub-CellList objects from an existing CellList. maybe subcells(celllist,[i,j,k])


# Project positions into a box assuming full periodicity
_getperiodicprojectedpos(pos,box) = mod.(pos[1], norm.(box))
# Get cell ids corresponding to positions
_getnormedcellids(pos,cshape,cnums,j) = 1 .+ Int.(mod.(ceil.(pos[j] ./ norm.(cshape)),cnums))

# Create a CellList our of an AtomsBase.FastSystem input
# Probably swap to kwargs once we know what we want to support
function CellList(data::AtomsBase.FastSystem, cutoff, cell_size_suggestion::AbstractVector{T}; boundaryconditions=:periodic, sortstyle = :unsorted) where T
    @assert all(in(cell_size_suggestion .≥ cutoff),true) "Cell size suggestions cannot be smaller than cutoff."
    @assert cutoff > zero(T) "Cutoff must be positive."
    box = data.bounding_box                                             # grab box dimensions
    d = size(box,1)                                                     # the dimension of the CellList is the same as that of the box
    cellNs = floor.(norm.(box) ./ cell_size_suggestion)                 # number of cells in each direction
    cellshape = box ./ cellNs                                           # compute cell axes & size since suggestion
                                                                        ## may not match integer multiple of box
    celllist = Cell{float(T)}.(collect(Iterators.product(range.(1,cellNs)...)),sortstyle) # generate array of empty & unsorted Cells
    coords = _getnormedcellids.((data.position,),(cellshape,),(cellNs,),1:length(data.position))    # Find out which cell ID each atom belongs in
    if boundaryconditions == :periodic
        # Account for periodicity in positions
        projectedpos = _getperiodicprojectedpos.(ustrip.(data.position),(ustrip.(box),)).*unit(eltype(box[1]))
    else
        error("Only :periodic boundary conditions are implemented so far.")
    end
    massivelyaddparticles!(celllist, projectedpos, data.atomic_symbol, coords)
    return CellList{typeof(celllist),typeof(cellshape),typeof(box),typeof(cutoff), size(box,1)}(celllist, cellshape, box, cutoff, boundaryconditions)
end

Base.show(io::IO, cl::CellList) = print(io, "CellList with ", ncells(cl), " cells")
Base.show(io::IO, ::MIME"text/plain", cl::CellList) = print(io, "CellList with ", size(cl), " cells in a " , cellgrid(cl), " grid.")

"""
size(cl::CellList)

The size of a CellList is the number of cells in the cell list.
"""
Base.size(cl::CellList) = prod(size(cl.celllist))

"""
natoms(cl)

The total number of atoms contained in a Cell, CellList or generic array of Cells.
"""
natoms(cl::CellList) = sum(length.(getfield.(cl.celllist, :position)))
natoms(cl::Cell) = length(getfield(cl, :position))
natoms(cl::Array{<:Cell}) = sum(length.(getfield.(cl, :position)))
"""
cellgrid(cl::CellList)

The grid layout of cells in a CellList or generic array of Cells.
"""
cellgrid(cl::CellList) = size(cl.celllist)

"""
massivelyaddparticles!(celllist::CellList, positions, elementtypes, cellcoords)

Add lots of atoms based on `positions` to `celllist` based on a list of cell IDs in `cellcoords`. `elementtypes` contains the element labels.
"""
function massivelyaddparticles!(celllist, positions, elementtypes, cellcoords)
    @inbounds for j = 1:length(positions)
        push!(celllist[cellcoords[j]...].position, positions[j]) # append atom
        push!(celllist[cellcoords[j]...].elementtypes, String.(elementtypes[j])) # append element of atom
        push!(celllist[cellcoords[j]...].originalindex, j) # append original index of atom in list
    end
    # adding new particles potentially broke sorting, so we have to resort the cells in the list
    cellsort!.(celllist)
    celllist
end

"""
cellsort!(cell::Cell, method)

Currently supported methods are :alphaelement (alphabetic sorting scheme by label names, e.g. atomic symbols) and :unsorted. Calling only `cellsort!(cell::Cell)`` will use the sorting scheme in `cell.sorttype`. Sorting a cell with a sorting scheme other than its own will update `cell.sorttype` to stay consistent with the new scheme.
"""
function cellsort!(cell::Cell, method::Symbol)
    cell.sorttype = method
    cellsort!(cell)
    cell
end
function cellsort!(cell::Cell)
    # TODO: I think this ends up with non-contiguous list in memory, possibly need a deepcopy here
    if cell.sorttype == :alphaelement
        perm = sortperm(cell.elementtypes)
        cell.position = cell.position[perm]
        cell.elementtypes = cell.elementtypes[perm]
        cell.originalindex = cell.originalindex[perm]
    elseif cell.sorttype == :unsorted
    else
        error("Sort style not implemented!")
    end
    # TODO: Implement more sorting styles.
    cell
end