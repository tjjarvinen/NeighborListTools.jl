

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

    tmp = sort_data(r, spc, cutoff)

    return tmp
end


form_cell_list(sys, cutoff; use_fp32=true) = form_cell_list(Array, sys, cutoff; use_fp32=use_fp32)


## 

@kernel function form_data(out, @Const(R), @Const(natom), @Const(cutoff), @Const(dcell))
    I =  @index(Global)
    @inbounds r =  R[I]


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


function sort_data(R, spc, cutoff, ncell)
    # sorting function
    f_sort(x,y) = x[2] < y[2]

    backend = get_backend(R)
    @assert backend == get_backend(spc)

    T = Tuple{Int32, Int32, eltype(spc), eltype(R)}
    data = similar(R, T)

    kernel = form_data(backend)
    kernel(data, R, spc, cutoff; ndrange=length(R))

    sort!(data; lt=f_sort)

    return data
end