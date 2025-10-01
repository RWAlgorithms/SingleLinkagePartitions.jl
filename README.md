# SingleLinkagePartitions
Minimalist single-linkage clustering.

# Documentation
The [documentation](https://rwalgorithms.github.io/SingleLinkagePartitions.jl/) for an overview and examples.

# Quick-start

```julia
import Random

import SingleLinkagePartitions as SL
using LinearAlgebra

T = Float64

#### singlet-linkage clustering.

# This is the set of points that we want to partition.
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

# # SLINK

# specification of point cloud, distance.
N = length(X)
X_mat = reshape(collect(Iterators.flatten(X)), length(X[begin]), length(X)) # alternatively, do X_mat = SL.vecs_to_mat(X)
dist_callable = SL.EuclideanDistance()

# Allocate state/buffer.
s = SL.SLINKState(T, N)

# run the SLINK algorithm. This computes the essential data required to construct the partition tree.
SL.slink!(s, dist_callable, X_mat)
# SL.slink!(s, dist_callable, X) # This would also work.

# construct the partition associated with level `num_fuses`. The possition levels take value from 0 to SL.get_cardinality(s), which is length(X) - 1.
num_fuses = 1
selected_partition = SL.construct_partition(s, num_fuses)

# construct all partitions, indexed by i := level - 1.
partition_tree = SL.construct_partition_tree(s)
```

The contents are:
```
julia> partition_tree = SL.construct_partition_tree(s)
10-element Memory{Memory{Memory{Int64}}}:
 [[1], [2], [3], [4], [5], [6], [7], [8], [9], [10]]
 [[1], [2], [3], [4], [5], [6], [8], [9], [7, 10]]
 [[2], [3], [4], [5], [6], [8], [1, 9], [7, 10]]
 [[2], [3], [4], [5], [8], [1, 9], [6, 7, 10]]
 [[2], [3], [5], [8], [1, 9], [4, 6, 7, 10]]
 [[3], [5], [8], [1, 2, 9], [4, 6, 7, 10]]
 [[5], [8], [1, 2, 3, 9], [4, 6, 7, 10]]
 [[5], [1, 2, 3, 8, 9], [4, 6, 7, 10]]
 [[1, 2, 3, 8, 9], [4, 5, 6, 7, 10]]
 [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]

julia> selected_partition
9-element Memory{Memory{Int64}}:
 [1]
 [2]
 [3]
 [4]
 [5]
 [6]
 [8]
 [9]
 [7, 10]
```

## Reduce points
We can use single-linkage partition tree to combine a set of points `X` in a part (e.g., cluster) of a partition. To select the partition given some suggested (this won't always be respected by the output) tolerance, `atol`, and zero tolerance for distance between duplicate points, `zero_tol`.

`level_trait` determines how we pick a level in the partition tree, which corresponds to picking a partition.

Each part in the selected partition, `partition_r`, generates a point in the output reduced point set `Xc`, and estimated variance `vs`. `center_trait` determines how we create a representative point given the points in a part, and we do this for each part in the partition to get the reduced point set.

```julia
# Let's purposely make the last point in the set to be very similar to the first point. We'll try to merge these two points in this example.
X[end] = X[1] .+ 0.00001 * randn(rng, 3)

atol = convert(T, 1.0)

# possible options for selecting a partition from the partition tree in `s`.
level_trait = SL.UseSLDistance(atol)
level_trait = SL.UseMaxDeviation(atol)
level_trait = SL.UseCumulativeSLDistance(atol)

# possible options for specifying the merged point for the points in a part of the partition.

# choose the closest input point in a part as the merged point for that part.
merge_trait = SL.UseProximitytoMean(copy_trait)

# choose the mean of a part as the merged point for that part.
merge_trait = SL.UseMean()

# choose the point that has the corresponding highest score array, `score`, as the merged point for the part it belongs to.
score = rand(rng, T, N)
merge_trait = SL.UseScore(copy_trait, SL.UseMaximum(), score)
# run `subtypes(SL.ScoreTrait)` to see a list of options that are similar to `SL.UseMaximum`.

# ## Reduce the input point set, `X`.
Xc0, vs0, partition_r0 = SL.reduce_pts(
    level_trait,
    merge_trait,
    dist_callable,
    X,
)

```

We could also pass in additional point sets to reduce, with the reduction based on the point distances in the primary point set, `X`.
```julia
# reduce both `X` and companion input `y`.
y = randn(rng, Complex{Float32}, N)
Xc1, vs_X1, yc, vs_y, partition_r1 = SL.reduce_pts(
    level_trait,
    merge_trait,
    dist_callable,
    X,
    y,
)

# reduce both `X` and companion input `y_set`.
y_set = collect(randn(rng, Complex{Float32}, N) for _ in 1:3)
Xc, vs_X, yc_set, vs_y_set, partition_r2 = SL.reduce_pts(
    level_trait,
    merge_trait,
    dist_callable,
    X,
    y_set,
)
```

Check if the reduced primary point set is the same as the one without companion inputs:
```julia
@assert norm(Xc0 - Xc1) < eps(T) * 100
@assert norm(Xc - Xc1) < eps(T) * 100

@assert norm(vs0 - vs_X1) < eps(T) * 100
@assert norm(vs_X - vs_X1) < eps(T) * 100
```

For each part in a given parition, we can also get the maximum magnitude deviation from the mean for each dimension.
```julia
ds_X = SL.computedeviationdims(X, partition_r)
```

### Routines based on `reduce_pts`
This package implements `avg_duplicates` and `replace_duplicates`, which are built on `reduce_pts`.

```julia
# average ponts that are very close together. `SL.UseSLDistance(a_tol)` is used to select a partition before averageing the points in each part.
X_avg, y_avg = SL.avg_duplicates(X, y, atol)
@assert abs((y[1] + y[end]) / 2 - y_avg[end]) < eps(T) * 100

# Similar to avg_duplicates, but does not use the mean as the representative merged point. Instead, use the point with the corresponding lowest score in y as the merged representative.
scores = randn(rng, T, length(X))
X_replaced, y_replaced = SL.replace_duplicates(copy_trait, SL.UseMaximum(), X, scores, atol)
```

# Citation
If you find this software useful, please cite this repository and code via the *Cite This Repository* button on the GitHub webpage of this repository.

# License
This project is licensed under the GNU Affero General Public License v3.0 license; see the `LICENSE` file for details. Individual source files may contain the following tag instead of the full license text:
```
SPDX-License-Identifier: AGPL-3.0-only
```

Using SPDX enables machine processing of license information based on the SPDX License Identifiers and makes it easier for developers to see at a glance which license they are dealing with.
