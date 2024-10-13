struct PairIterator{TCL, TI}
    clist::TCL
    indices::TI
    npairs::Int
end

Base.IteratorSize(::PairIterator) = Base.SizeUnknown()


function give_site_iterators(cl::CellList, n::Int)
    return map( chunks( CartesianIndices(cl.cell_indices); n=n)  ) do c
        return CellIterator(cl, c, true)
    end
end


function give_pair_iterators(cl::CellList, n::Int; npairs::Int=0)
    tmp =  map( chunks( CartesianIndices(cl.cell_indices); n=n) ) do c
        CellIterator(cl, c, npairs)
    end
    return tmp
end


function _add_lists(a,b)
    return (
        r = vcat(a.r, b.r),
        origin_indx1 = vcat(a.origin_indx1, b.origin_indx1),
        origin_indx2 = vcat(a.origin_indx2, b.origin_indx2),
        species1 = vcat(a.species1, b.species1),
        species2 = vcat(a.species2, b.species2),
    )
end

function _cut_last_ones(buffer, i)
    return (
        r = buffer.r[i:end],
        origin_indx1 = buffer.origin_indx1[i:end],
        origin_indx2 = buffer.origin_indx2[i:end],
        species1 = buffer.species1[i:end],
        species2 = buffer.species2[i:end],
    )
end


function  Base.iterate(pi::PairIterator)
    new = iterate(pi.indices)
    isnothing(new) && return nothing

    tmp_item, istate = new

    tmp = get_pairlist(pi.clist, tmp_item)
    if length(tmp.r) < pi.npairs
        new = iterate(pi.indices, istate)
        while length(tmp.r) <= pi.npairs && ! isnothing(new)
            tmp_item, istate = new
            t = get_pairlist(pi.clist, tmp_item)
            tmp = _add_lists(tmp, t)
            new = iterate(pi.indices, istate)
        end
        if length(tmp.r) < pi.npairs  
            # not enough pairs for needed ammount
            return tmp, (istate=nothing, buffer=nothing, i=nothing)
        end
    end
    state = (istate=istate, buffer=tmp, i=pi.npairs)
    item = (
        r = tmp.r[1:pi.npairs],
        origin_indx1 = tmp.origin_indx1[1:pi.npairs],
        origin_indx2 = tmp.origin_indx2[1:pi.npairs],
        species1 = tmp.species1[1:pi.npairs],
        species2 = tmp.species2[1:pi.npairs],
    )
    return item, state
end

function  Base.iterate(pi::PairIterator, state)
    isnothing(state.istate) && return nothing

    if length(state.buffer.r)-state.i > pi.npairs
        # We have enough buffer
        i = state.i + pi.npairs
        new_state = (istate=state.istate, buffer=state.buffer, i=i)
        buffer = state.buffer
        indx = (state.i+1):i 
        item = (
            r = buffer.r[indx],
            origin_indx1 = buffer.origin_indx1[indx],
            origin_indx2 = buffer.origin_indx2[indx],
            species1 = buffer.species1[indx],
            species2 = buffer.species2[indx],
        )
        return item, new_state
    end

    new = iterate(pi.indices, state.istate)
    if isnothing(new)
        if length(state.buffer.r)-state.i > 0
            # return what we have
            i = state.i + 1
            new_state = (istate=nothing, buffer=nothing, i=nothing)
            buffer = state.buffer
            item = (
                r = buffer.r[i:end],
                origin_indx1 = buffer.origin_indx1[i:end],
                origin_indx2 = buffer.origin_indx2[i:end],
                species1 = buffer.species1[i:end],
                species2 = buffer.species2[i:end],
            )
            return item, new_state
        end
        return nothing
    end

    tmp_item, istate = new
    tmp = get_pairlist(pi.clist, tmp_item)
    tmp_buff = _cut_last_ones(state.buffer, state.i+1)
    tmp = _add_lists(tmp_buff, tmp)

    if length(tmp.r) < pi.npairs
        new = iterate(pi.indices, istate)
        while length(tmp.r) <= pi.npairs && ! isnothing(new)
            tmp_item, istate = new
            t = get_pairlist(pi.clist, tmp_item)
            tmp = _add_lists(tmp, t)
            new = iterate(pi.indices, istate)
        end
        if length(tmp.r) < pi.npairs  
            # not enough pairs return what we have
            return tmp, (istate=nothing, buffer=nothing, i=nothing)
        end
    end
    state = (istate=istate, buffer=tmp, i=pi.npairs)
    item = (
        r = tmp.r[1:pi.npairs],
        origin_indx1 = tmp.origin_indx1[1:pi.npairs],
        origin_indx2 = tmp.origin_indx2[1:pi.npairs],
        species1 = tmp.species1[1:pi.npairs],
        species2 = tmp.species2[1:pi.npairs],
    )
    return item, state
end

##


struct CellIterator{TCL, TI}
    clist::TCL
    indices::TI
    sites::Bool
end

Base.length(cl::CellIterator) = length(cl.indices)


function Base.iterate(cl::CellIterator)
    new = iterate(cl.indices)
    isnothing(new) && return nothing
    item, state = new
    if cl.sites
        tmp = get_site_pairlist(cl.clist, item)
    else
        tmp = get_pairlist(cl.clist, item)
    end
    return tmp, state
end


function Base.iterate(cl::CellIterator, state)
    new = iterate(cl.indices, state)
    isnothing(new) && return nothing
    item, state = new
    if cl.sites
        tmp = get_site_pairlist(cl.clist, item)
    else
        tmp = get_pairlist(cl.clist, item)
    end
    return tmp, state
end