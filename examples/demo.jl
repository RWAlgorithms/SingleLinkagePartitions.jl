import Random
rng = Random.Xoshiro(0)

# import SingleLinkagePartitions
# SL = SingleLinkagePartitions


#### singlet-linkage clustering.
T = Float64
copy_trait = SL.MakeCopy()

# # replacement by extrema score.
U0 = collect(ones(1) for _ in 1:3)
U0[end][1] = 3
y0 = [1.1; 1.2; 1.3]
Uc, yc = SL.replace_duplicates(copy_trait, SL.UseMinimum(), U0, y0, 1.0e-5)
println("replace_duplicates() with UseMinimum")
@show U0, y0
@show Uc, yc
println()

U0 = collect(ones(1) for _ in 1:3)
U0[end][1] = 3
y0 = [1.1; 1.2; 1.3]
Uc, yc = SL.replace_duplicates(copy_trait, SL.UseMaximum(), U0, y0, 1.0e-5)
println("replace_duplicates() with UseMaximum")
@show U0, y0
@show Uc, yc
println()


# # from the README.md.
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
D = 3

# Let's purposely make the last point in the set to be very similar to the first point. We'll try to merge these two points in this example.
X[end] = X[1] .+ 0.00001 * randn(rng, 3)

# # SLINK

# specification of point cloud, distance.
N = length(X)
X_mat = reshape(collect(Iterators.flatten(X)), length(X[begin]), length(X)) # alternatively, do X_mat = SL.vecs_to_mat(X)
dist_callable = SL.EuclideanDistance()

# Allocate state/buffer.
s = SL.SLINKState(T, N)

# run the SLINK algorithm. This computes the essential data required to construct the partition tree.
SL.slink!(s, dist_callable, X_mat)

# construct the partition associated with level `num_fuses`. The possition levels take value from 0 to SL.get_cardinality(s), which is length(X) - 1.
num_fuses = 1
selected_partition = SL.construct_partition(s, num_fuses)

# construct all partitions, indexed by i := level - 1.
partition_tree = SL.construct_partition_tree(s)

# # Reduce points
# Select a partition from the SL tree based on level_trait <: SL.LevelOption.
# For each part in the partition, replace all the points in that part with a representative point. The user use`merge_trait <: SL.PartRepOption` to specify the strategy for selecting the representative point.

atol = convert(T, 1.0e-4)

# possible options for selecting a partition from the partition tree in `s`.
level_trait = SL.UseSLDistance(atol)
level_trait = SL.UseMaxDeviation(atol)
level_trait = SL.UseCumulativeSLDistance(atol)

# possible options for specifying the merged point for the points in a part of the partition.

# choose the closest input point in a part as the merged point for that part.
merge_trait = SL.UseProximitytoMean(copy_trait)

# choose the mean of a part as the merged point for that part.
merge_trait = SL.UseMean()

# choose the point that has the corresponding highest score array, denoted `score` in the following code, as the merged point for the part it belongs to.
score = rand(rng, T, N)
merge_trait = SL.UseScore(copy_trait, SL.UseMaximum(), score)

# ## Reduce the input point set, `X`.
Xc0, vs0, partition_r0 = SL.reduce_pts(
    level_trait,
    merge_trait,
    dist_callable,
    X,
)

println("The last point of Y is to be merged because the last slot of partition_r0 has a non-singleton part.")
@show partition_r0

# reduce both `X` and companion point set, `y`. This version is when the points in `y` is 1-D.
y = randn(rng, Complex{Float32}, N)
Xc1, vs_X1, yc, vs_y, partition_r2 = SL.reduce_pts(
    level_trait,
    merge_trait,
    dist_callable,
    X,
    y,
)

# Companion point set is multi-dim.
D_y = 2
Y = collect(randn(rng, Complex{Float32}, D_y) for _ in 1:N)
Xc3, vs_X3, Yc3, vs_Y3, partition_r3 = SL.reduce_pts(
    level_trait,
    merge_trait,
    dist_callable,
    X,
    Y,
)

# make sure the reduced primary input set is the same, between the companion input versions of `reduce_pts`.
@assert all(partition_r0 .== partition_r2)
@assert all(partition_r3 .== partition_r2)

@assert norm(Xc0 - Xc1) < eps(T) * 100
@assert norm(Xc0 - Xc3) < eps(T) * 100
@assert norm(vs0 - vs_X3) < eps(T) * 100

# ## Routines built on reduce_pts

# average ponts that are very close together. `SL.UseSLDistance(a_tol)` is used to select a partition before averageing the points in each part.
X_avg, y_avg = SL.avg_duplicates(X, y, atol)

# Similar to avg_duplicates, but does not use the mean as the representative merged point. Instead, use the point with the corresponding lowest score in y as the merged representative.
scores = randn(rng, T, length(X))
X_replaced, y_replaced = SL.replace_duplicates(copy_trait, SL.UseMaximum(), X, scores, atol)

nothing
