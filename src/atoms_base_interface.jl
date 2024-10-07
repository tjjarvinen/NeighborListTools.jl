

function CellList(sys, cutoff::Unitful.Length)
    form_cell_list(sys, cutoff)
end


function form_cell_list(TA, sys, cutoff; use_fp32=true)
    box = bounding_box(sys)
    # only orthorhombic cell allowed
    @argcheck box[1][2] ≈ box[1][3] ≈ zero(box[1][3])
    @argcheck box[2][1] ≈ box[2][3] ≈ zero(box[2][3])
    @argcheck box[3][1] ≈ box[3][2] ≈ zero(box[3][1])

    r = position(sys, :)
    if use_fp32
        r = map(r) do rᵢ
            Float32.(rᵢ)
        end
    end
    r = TA(r)
    spc = TA(species(sys, :))
    pbc = periodicity(sys)

    indx = similar(spc, Int32)
    cell_indx = similar(indx)

    cell_diag = SVector( box[1][1], box[2][2], box[3][3] )
    if use_fp32
        cell_diag = Float32.(cell_diag)
    end

    # Fit cutoff to cell boundaries
    ncells = floor.(Int32, cell_diag / cutoff)
    coff = cell_diag ./ ncells

    tmp = sort_data(r, spc, coff, cell_diag)

    cells = move_data_to_cells!(indx, cell_indx, r, spc, tmp, ncells)

    return CellList(
        box,
        cells.cells,
        cutoff,
        cells.index,
        cells.pos,
        cells.species
    )
end


form_cell_list(sys, cutoff; use_fp32=true) = form_cell_list(Array, sys, cutoff; use_fp32=use_fp32)


## 

@kernel function form_data(out, @Const(R), @Const(natom), @Const(cutoff), @Const(dcell))
    I =  @index(Global)
    @inbounds r =  R[I]

    # wrap atoms inside the box
    @inbounds r1 = r[1] - (div(r[1], dcell[1]) -1) * dcell[1]
    @inbounds r1 = r1 - div(r1, dcell[1]) * dcell[1]
    @inbounds r2 = r[2] - (div(r[2], dcell[2]) -1) * dcell[2]
    @inbounds r2 = r2 - div(r2, dcell[2]) * dcell[2]
    @inbounds r3 = r[3] - (div(r[3], dcell[3]) -1) * dcell[3]
    @inbounds r3 = r3 - div(r3, dcell[3]) * dcell[3]

    r = SVector(r1, r2, r3)

    # Calculate cell index
    @inbounds tia = floor(r[1]/cutoff[1])
    ia = unsafe_trunc(Int32, tia)
    @inbounds tib = floor(r[2]/cutoff[2])
    ib = unsafe_trunc(Int32, tib)
    @inbounds tic = floor(r[3]/cutoff[3])
    ic = unsafe_trunc(Int32, tic)

    iabc = (ia+1) + 1_000*(ib+1) + 1_000_000*(ic+1)

    @inbounds out[I] = (I, iabc, natom[I], r)
end


function sort_data(R, spc, cutoff, dcell)
    # sorting function
    f_sort(x,y) = x[2] < y[2]

    backend = get_backend(R)
    @assert backend == get_backend(spc)

    T = Tuple{Int32, Int32, eltype(spc), eltype(R)}
    data = similar(R, T)

    kernel = form_data(backend)
    kernel(data, R, spc, cutoff, dcell; ndrange=length(R))

    sort!(data; lt=f_sort)

    return data
end


function move_data_to_cells!(indx, cell_indx, r, spc,  data, ncells)
    @argcheck length(r) == length(spc) == length(data)

    # This needs to be optimized
    map!( x->x[1], indx, data )
    map!( x->x[2], cell_indx, data )
    map!( x->x[3], spc, data )
    map!( x->x[4], r, data )

    cells = similar(r, UnitRange{Int}, (ncells...) )

    for ix in axes(cells,1), iy in axes(cells, 2), iz in axes(cells, 3)
        cind = ix + iy*1000 + iz*1_000_000

        # NOTE this does not work with GPU yet
        # https://github.com/anicusan/AcceleratedKernels.jl has it
        # but has not yet been registered yet
        cells[ix,iy,iz] = searchsorted(cell_indx, cind)
    end

    return (pos=r, species=spc, index=indx, cell_index=cell_indx, cells=cells)
end
