# NeighborListTools

[![Build Status](https://github.com/tjjarvinen/NeighborListTools.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tjjarvinen/NeighborListTools.jl/actions/workflows/CI.yml?query=branch%3Amain)

This is an experimental neighbor list generating package.
That is designed to also work with GPU arrays.

Currently usable for cases with orthogonal cell vectors and cell larger than cutoff radius.

## Example

```julia
using AtomsBuilder
using NeighborListTools
using Unitful

cutoff = 9.0u"Å"
sys = bulk(:Ar, cubic=true) * 15

cl = CellList(sys, cutoff)

# Cell list size
size(cl)

# Get cell (1,1,1)
cl[1,1,1]

# Get pairlist for cell (1,1,1)
# this is pair interaction list
get_pairlist(cl, 1,1,1)

# this is site potential list
pl = get_site_pairlist(cl, 1,1,1)

# distance vector for the pairs
pl.r

# Vectors of species of particles
pl.species1
pl.species2

# Vectors of indices in the original data
pl.origin_indx1
pl.origin_indx2

# Vectors of indices in the cell list
pl.indx1
pl.indx2

# Go through the whole pair list
for i in eachindex(cl)
    pl = get_pairlist(cl, i)
    # do something
end
```

## Iterator interface

Note, this is experimental.

```julia

# Give 4 pair interaction iterators
iterators = give_pair_iterators(cl, 4)
# or give 4 site potential pair iterators
iterators = give_site_iterators(cl, 4)

Threads.@threads for iter in iterators
    for neighbors in iter
        # do something
        nothing
    end
end

```



## GPU support

GPU support is not yet fully working.

You can create GPU CellList by

```julia
using AtomsBuilder
using NeighborListTools
using Unitful
using CUDA

cutoff = 9.0u"Å"
sys = bulk(:Ar, cubic=true) * 15

cl = CellList(CuArray, sys, cutoff)

pl = get_pairlist(cl, 1,1,1)
```

Iterator interface does not work with GPU.
GPU is also very slow compared to CPU due to kernel lauch latency.


## SIMD

Output now includes vectors that are build with [SIMD](https://github.com/eschnett/SIMD.jl) use in mind.

```julia
using AtomsBuilder
using NeighborListTools
using Unitful
using SIMD

cutoff = 9.0u"Å"
sys = bulk(:Ar, cubic=true) * 15

cl = CellList( sys, cutoff)

pl = get_pairlist(cl, 1,1,1)

# x-,y- and z-coordinates for the particles
# these are intended to be used with vload from SIMD.jl
pl.x
pl.y
pl.z

# this mask for filtering out at the end of x-, y- and z-arrays
# that are needed to complete the vector length
pl.mask


# Lets have Lennard-Jones energy
function my_energy(x,y,z, c12=Float32(0.3), c6=Float32(2.0))
    r2 = x^2 + y^2 + z^2
    r6 = r2^3
    r12 = r6^2
    return c12/r12 - c6/r6
end


E_tot = Vec{16,Float32}(0)
for i in 1:16:length(pl.x)
    x = vload(Vec{16, Float32}, pl.x, i)
    y = vload(Vec{16, Float32}, pl.y, i)
    z = vload(Vec{16, Float32}, pl.z, i)
    mask = vload(Vec{16, Float32}, pl.mask, i)
    e = my_energy(x,y,z)  # this is now vectorized
    E_tot = muladd(mask, e, E_tot) # this is also vectorized
end
sum(E_tot) # total energy
```