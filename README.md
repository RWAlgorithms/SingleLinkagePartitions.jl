# SingleLinkagePartitions
Minimalist single-linkage clustering. Follows the algorithm by Sibson.
Thanks to Dr. Veronika Strnadova-Neeley for her online notes on this topic.

# Quick-start

```julia
import Random
Random.seed!(25)

import SingleLinkagePartitions
using LinearAlgebra

#### singlet-linkage clustering.

# generate data randomly:
# D = 3
# N = 10
# X = collect( randn(D) for _ = 1:N )

# or, use this default set of points.
X = [
    [-0.40243248293794137, 0.8540414903329187, -0.6651248667822778],
    [-0.9754736537655263, -0.341379056315598, -1.0410555755312705],
    [-1.0496381529964869, -1.0289732912432703, -0.4305269991795779],
    [0.7262044632757468, 0.8138894370909177, 0.6104189261116074],
    [2.0501294946950264, 0.18095967976913974, 0.9406747232875855],
    [1.0407043988018494, -0.14776493165237461, -0.8737149501414327],
    [1.0484097740458416, 0.7379871044247477, -0.02494318852134621],
    [-0.32639477363891256, -1.45405586112584, 0.5104603413606108],
    [-0.6283556049853254, 0.35921840490046464, -1.1166717373759707],
    [1.2421363428579315, 0.47437350434528236, -0.5869506255304089],
]


# specify metric function.
metricfunc = (xx,yy)->norm(xx-yy)

# get single linkage partitions. Pass in `early_stop_distance = Inf` to generate all partitions until we have the singleton partition.
distance_set, partition_set = SingleLinkagePartitions.runsinglelinkage(
    X,
    metricfunc;
    early_stop_distance = Inf,
)

Y, status_flag = SingleLinkagePartitions.mergepoints(
    X,
    metricfunc;
    tol = 0.6,
)


```

The output:
```
julia> distance_set
10-element Vector{Float64}:
 0.0
 0.6502879922496319
 0.7069552456099949
 0.7140482053065045
 0.7164233686190574
 0.7855224671411315
 0.9225136028232725
 1.2606479889562474
 1.4987127483818177
 1.590045425866602

julia> partition_set
10-element Vector{Vector{Vector{Int64}}}:
 [[1], [2], [3], [4], [5], [6], [7], [8], [9], [10]]
 [[1], [2], [3], [4], [5], [6], [8], [9], [10, 7]]
 [[2], [3], [4], [5], [6], [8], [9, 1], [10, 7]]
 [[2], [3], [4], [5], [8], [9, 1], [10, 7, 6]]
 [[2], [3], [5], [8], [9, 1], [10, 7, 6, 4]]
 [[3], [5], [8], [9, 1, 2], [10, 7, 6, 4]]
 [[5], [8], [9, 1, 2, 3], [10, 7, 6, 4]]
 [[5], [9, 1, 2, 3, 8], [10, 7, 6, 4]]
 [[9, 1, 2, 3, 8], [10, 7, 6, 4, 5]]
 [[10, 7, 6, 4, 5, 9, 1, 2, 3, 8]]

 julia> Y
10-element Vector{Vector{Float64}}:
 [-0.40243248293794137, 0.8540414903329187, -0.6651248667822778]
 [-0.9754736537655263, -0.341379056315598, -1.0410555755312705]
 [-1.0496381529964869, -1.0289732912432703, -0.4305269991795779]
 [0.7262044632757468, 0.8138894370909177, 0.6104189261116074]
 [2.0501294946950264, 0.18095967976913974, 0.9406747232875855]
 [1.0407043988018494, -0.14776493165237461, -0.8737149501414327]
 [1.0484097740458416, 0.7379871044247477, -0.02494318852134621]
 [-0.32639477363891256, -1.45405586112584, 0.5104603413606108]
 [-0.6283556049853254, 0.35921840490046464, -1.1166717373759707]
 [1.2421363428579315, 0.47437350434528236, -0.5869506255304089]

julia> status_flag
true
```

If `mergepoints()` returns a `false` value for `status_flag`, then an unexpected error occured, and the pair-wise distances (as computed by `metricfunc`) returned merged points `Y` are not all smaller than `tol`.