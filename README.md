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

# distance vector for pairs
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
```

Trying to calculate pairlist fails though.
