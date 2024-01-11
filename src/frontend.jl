
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