# SingleLinkagePartitions
Minimalist single-linkage clustering. 

# Documentation
The [documentation](https://rwalgorithms.github.io/SingleLinkagePartitions.jl/) for an overview and examples.

# Quick-start

```julia
import Random
Random.seed!(25)

import SingleLinkagePartitions as SL
using LinearAlgebra

T = Float64

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
N = length(X)

metric = SL.geteuclideanmetric()
pt = SL.computesl(metric, X)

distances_set = SL.getdistances(pt)
partition_set = SL.generateallpartitions(pt)
```

The output:
```
julia> distance_set = SL.getdistances(pt)
9-element Vector{Float64}:
 0.6502879922496316
 0.7069552456099949
 0.7140482053065044
 0.7164233686190572
 0.7855224671411318
 0.9225136028232729
 1.2606479889562476
 1.4987127483818177
 1.5900454258666015

julia> partition_set = SL.generateallpartitions(pt)
10-element Vector{Vector{Vector{Int64}}}:
 [[1], [2], [3], [4], [5], [6], [7], [8], [9], [10]]
 [[1], [2], [3], [4], [5], [6], [7, 10], [8], [9]]
 [[1, 9], [2], [3], [4], [5], [6], [7, 10], [8]]
 [[1, 9], [2], [3], [4], [5], [6, 7, 10], [8]]
 [[1, 9], [2], [3], [4, 6, 7, 10], [5], [8]]
 [[1, 2, 9], [3], [4, 6, 7, 10], [5], [8]]
 [[1, 2, 3, 9], [4, 6, 7, 10], [5], [8]]
 [[1, 2, 3, 8, 9], [4, 6, 7, 10], [5]]
 [[1, 2, 3, 8, 9], [4, 5, 6, 7, 10]]
 [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
```


## Reduce points
We can use single-linkage partition tree to combine a set of points `X` in a part (e.g., cluster) of a partition. To select the partition given some suggested (this won't always be respected by the output) tolerance, `atol`, and zero tolerance for distance between duplicate points, `zero_tol`.

`level_trait` determines how we pick a level in the partition tree, which corresponds to picking a partition.

Each part in the selected partition, `partition_r`, generates a point in the output reduced point set `Xc`, and estimated variance `vs`. `center_trait` determines how we create a representative point given the points in a part, and we do this for each part in the partition to get the reduced point set.

```julia
atol = convert(T, 1.0)
level_trait = SL.UseSLDistance(atol) # `atol` is compared aaginst`distance_set to pick a level.

# other trais.
# level_trait = SL.UseCumulativeSLDistance() # atol is compared against cumsum(distance_set) to pick a level.

# For each part, pick the point that has the lowest corresponding score among the points in the part as the representative for that part.
score = randn(T, N)
center_trait = SL.UseScore(SL.UseMinimum(), score) # use SL.UseMaximum() to select based on the highest score.

# For each part, pick the point that has the closest Euclidean distance to the mean among the points in the part as the representative for that part. The mean is taken over all the points of the part.
# center_trait = SL.UseProximitytoMean() # when Xc must be a subset of X

# For each part, assign the mean over all the points in the part as the representative. This implies we are not using a point in the original point set `X` as the representative points, so `Xc` is not a subset of `X`.
# center_trait = SL.UseMean() # when Xc can be any pt from which X is constructed 

Xc, vs, partition_r = SL.reducepts(
    level_trait,
    center_trait,
    metric,
    X,
)
```

We could also pass in an accompanying scalar value for each point in `X`. These scalar values make up `y`, and can be a complex-valued or real-valued floating-point number that isn't of the same data type as the points in `X`. Whatever scheme we used to pick a representative for each part for `Xc`, we do for `yc`. Similarly, we provide a variance `vs_y` for `yc`, and a set of variances `vs_X` for `Xc`.
```julia

y = randn(Complex{Float32}, N)
Xc, vs_X, yc, vs_y, partition_r = SL.reducepts(
    level_trait,
    center_trait,
    metric,
    X,
    y,
)

```

For each part in a given parition, we can also get the maximum magnitude deviation from the mean for each dimension.
```julia
ds_X = SL.computedeviationdims(X, partition_r)
```

# Citation
You can use the *Cite This Repository* button below the *About* section on the GitHub webpage of this repository.
