
Random.seed!(25)

using Printf

PLT.close("all")
fig_num = 1

# # Set up data
T = Float64
D = 2;

# oracle set to one branch of Fermat's spiral
f = tt->(sqrt(tt) .* [cos(tt); sin(tt)]);

# sample data from the oracle to get our 2-D point set, X.
N = 10
ts = LinRange(0, 2*π, N)
X = collect( f(t) for t in ts ); # this is our data, the input point set for this demo.

# visualize the data.
plot_x = collect( X[n][begin] for n in eachindex(X) )
plot_y = collect( X[n][begin+1] for n in eachindex(X) )

PLT.figure(fig_num)
fig_num += 1
PLT.scatter(plot_x, plot_y, marker = "x")
PLT.axis("scaled")
PLT.title("Input point set")
PLT.gcf()

# ## Single-linkage clustering
# Single-linkage clustering is a method to generate a partition tree, which is a set of nested partitions.
# First, specify the distance metric.
metric = SL.geteuclideanmetric()

# ### Generate partition tree
pt = SL.computesl(metric, X)

# The minimum distances between the nested partitions.
distance_set = SL.getdistances(pt) # length(distance_set) == length(X) - 1.

# The nested partitions in the partition tree.
partition_set = SL.generateallpartitions(pt) # length(partition_set) == length(X).


# We can see that only a single part gets larger and larger. This is because the minimum distance between all parts for a partition at a given level is one that involves the large part. This is a known characteristic with single-linkage clustering, which may be something that is desirable or undesirable, depending on the application.

# ## Select a level (i.e. partition) via conditions
# The conditions implemented are: `UseMaxDeviation`, `UseSLDistance`, and `UseCumulativeSLDistance`.

# ### `UseMaxDeviation`: Deivation from part/cluster mean/centroid.
# See the terminology section.

# specify an allowed deviation distance.
max_dev = pt.w[end]/2 # we arbitrarily choose half of the largest single-linkage distance from the tree.
level_config = SL.UseMaxDeviation(max_dev)

level = SL.picklevel(level_config, pt, X) # pick level via binary search.
partition = SL.getpartition(pt, level)

# We see that the maximum deviation sequence when ordered by the nesting order of the partition tree (i.e. from level = 0 the leaf to level == getmaxlevel(pt) the root) is approximately monotonic.
max_ds = collect(
    SL.computemaxdeviation(X, pt, level)
    for level in 0:(SL.getmaxlevel(pt))
)
PLT.figure(fig_num)
fig_num += 1
PLT.plot(0:(SL.getmaxlevel(pt)), max_ds, "o")
#PLT.axis("scaled")
PLT.xlabel("Level index, in the nested order")
PLT.ylabel("The Maximum deviation of a partition")
PLT.title("Maximum deviation of all partitions in the tree")
PLT.gcf()

# Since it is approximately monotonic, we use a bracketed binary search to select a partition such that its maximum deviation is approximate the best match from all the partitions to the specified value `max_dev`.
# One can alternatively do a linear search over `max_ds` to pick a level, which is easy to implement yourself given `max_ds``.
level = SL.picklevel(level_config, pt, X) # pick level via binary search.
partition = SL.getpartition(pt, level) # instantiate the partition given the level and tree.

println("Target maximum deviation: ", max_dev)
println("Selected partition's maximum deviation: ", SL.computemaxdeviation(X, pt, level));

# visualize
marker_size = 100.5
tmp = @sprintf("%02i", max_dev)
title_string = "Partition selected via UseMaxDeviation"
ARGS = (partition, title_string, X, fig_num, marker_size) # hide
include("./snippet/plot_partition.jl")

# # Iterated single-linkage to reduce chaining
# If it is desirable to reduce chaining, then we can use an iterated version of picklevel(). So far, only `UseMaxDeviation` is implemented for iteration.

# We need to specify the iteration, namely a discount factor-like parameter that characterizes our acceptance of a selected partition in every iteration.
acceptance_factor = 0.99

# `acceptance_factor` is a value in (0,1). closest to 1 means almost only the part with the maximum dev is added to P.
# A larger value for `acceptance_factor`` tend to yield fewer parts in the final P.
# We now run the iteration
P, iters_ran = SL.iteratedsl(
    level_config,
    SL.geteuclideanmetric(),
    X;
    acceptance_factor = acceptance_factor,
    max_iter = 100,
)

max_dev_P = SL.computemaxdeviation(X, partition)

# visualize
partition = P
tmp = @sprintf("%02i", max_dev)
title_string = "Partition selected via iteratedsl"
ARGS = (partition, title_string, X, fig_num, marker_size) # hide
include("./snippet/plot_partition.jl")


nothing