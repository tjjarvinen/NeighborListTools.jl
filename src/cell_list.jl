
## Cell and all its methods

mutable struct Cell{TA, TI}
    "positions of atoms"
    positions::TA
    "indexes in original input"
    indexes::TI
    "number of atoms in the cell"
    natoms::Int   # natoms <= lenght(positions) - use to make view on the cell content
    function Cell(allog_length::Int; array_type=Array, element_type=SVector{3, Float64})
        pos = array_type{element_type}(undef, allog_length)
        ind = array_type{Int}(undef, allog_length)
        new{typeof(pos), typeof(ind)}(pos, ind, 0)
    end
end

Base.show(io::IO, c::Cell) = print(io, "Cell with ", length(c), " atoms")
Base.show(io::IO, ::MIME"text/plain", c::Cell) = print(io, "Cell with ", length(c), " atoms")

Base.length(c::Cell) = c.natoms


## CellList and its methods

mutable struct CellList{T, TC, N, ND} # N is number of different types; ND is number of dimensions
    cell_list::T
    cutoff::TC
    information::Dict{Symbol, Any}
    function CellList(base_cell::Cell, cutoff, cell_shape...; kwargs...)
        tmp_list = fill(base_cell, cell_shape) # Note these are copies of the original
        list = map( x->deepcopy(x), tmp_list)  # Now they are unique
        info = Dict(kwargs...)
        new{typeof(list), typeof(cutoff), cell_shape[1], length(cell_shape)-1}(list, cutoff, info)
    end
end

Base.show(io::IO, cl::CellList) = print(io, "CellList with ", size(cl), " cells")
Base.show(io::IO, ::MIME"text/plain", cl::CellList) = print(io, "CellList with ", size(cl), " cells")

Base.size(cl::CellList) = size(cl.cell_list)[2:end]

number_of_types(cl::CellList) = size(cl.cell_list, 1)