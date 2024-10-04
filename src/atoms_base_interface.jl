

function form_cell_list(TA, sys, cutoff; use_fp32=true)
    r = position(sys, :)
    if use_fp32
        r = map(r) do rᵢ
            Float32.(rᵢ)
        end
    end
    spc = species(sys, :)
    box = bounding_box(sys)
    pbc = periodicity(sys)
    cell_diag = SVector( box[1][1], box[2][2], box[3][3] )
    if use_fp32
        cell_diag = Float32.(cell_diag)
    end

    tmp = sort_data(r, spc, cutoff, cell_diag)

    return tmp
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
    @inbounds tia = ceil(r[1]/cutoff)
    ia = unsafe_trunc(Int32, tia)
    @inbounds tib = ceil(r[2]/cutoff)
    ib = unsafe_trunc(Int32, tib)
    @inbounds tic = ceil(r[3]/cutoff)
    ic = unsafe_trunc(Int32, tic)

    iabc = ic + 1_000*ib + 1_000_000*ia

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