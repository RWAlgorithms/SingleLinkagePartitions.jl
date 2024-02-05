
"""
```
getmaxlevel(pt::PartitionTree)::Int
```
Returns the maximum level of the partition tree.

A valid level for `pt` can take integer values in [0, getmaxlevel(pt)].
A level indexes a partition in `pt`.
"""
function getmaxlevel(pt::PartitionTree)::Int
    return getNedges(pt)
end

"""
    avgduplicates(X::Vector, y::Vector, zero_tol)

For any locations in X that are nearly identical (with respect to `zero_tol`), reduce the corresponding `y` entries by `mean`.
`zero_tol` should be set to a very small number, such as a multiple of `eps`.

Returns
- `Xc`: a reduced version of `X`.
- `yc`: a reduced version of `y`.
"""
function avgduplicates(X::Vector, y::Vector, zero_tol)

    # when Xc can be any pt from which X is constructed from.
    center_trait = UseMean() # reduce pts via mean().
    level_trait = UseSLDistance(zero_tol)

    Xc, _, yc, _ = reducepts(
        level_trait,
        center_trait,
        geteuclideanmetric(),
        X,
        y,
    )

    return Xc, yc
end