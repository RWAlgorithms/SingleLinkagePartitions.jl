using Random
Random.seed!(25)


T = Float64
Δc, xs, Fs = deserialize("data/rg_Float64")
@show length(Δc)

# # set up point set and max deviation from mean condition.
# This condition is such that each part in the selected partition has a maximum magnitude of max_dev from the mean of that part.

X = copy(Δc)

max_dev = 0.2
level_config = SL.UseMaxDeviation(max_dev)

metric = SL.geteuclideanmetric()

###### option 1: single run of SL to get partition tree, then binary search the tree to get the partition that satisfies the max deviation from mean condition.

pt = SL.computesl(metric, X) # get partition tree.
level = SL.picklevel(level_config, pt, X) # pick level via binary search.
partition = SL.getpartition(pt, level) # instantiate the partition given the level and tree.
@show length(partition)

# compute the maximum magnitude deviateion from the mean.
ds_X0 = SL.computedeviationdims(X, partition)
max_ds0 = collect(maximum(ds_X0[k]) for k in eachindex(ds_X0))
@show maximum(max_ds0)
println()

##### option 2: repeated runs of SL, each computing a partition tree.
# At the end of each run, move the parts from the current partition tree to the output partition, P, if it has within acceptance_factor*Z amount of magnitude deviation from the mean,
# where Z is the maximum magnitude deviation out of all the parts in the current partition tree.
# The points not selected in the run is assembled into a new input X, and repeated.
# This discourages "chaining", so the deviations within parts tend to be more uniform, while option 1 would have chaining, where the larger parts have a large deviation from its mean when compared to the smaller ones.

acceptance_factor = 0.99 # a value in (0,1). closest to 1 means almost only the part with the maximum dev is added to P.
# A larger value tend to yield fewer parts in the final P.

P, iters_ran = SL.iteratedsl(
    level_config,
    SL.geteuclideanmetric(),
    X;
    acceptance_factor = acceptance_factor,
    max_iter = 100, # how many iterations we allow.
)
@show length(P)


ds_X = SL.computedeviationdims(X, P)
max_ds = collect(maximum(ds_X[k]) for k in eachindex(ds_X))

@show maximum(max_ds)
println()

# sanity check: P should contain unique entries that take value from 1:length(X).
P_flat = collect(Iterators.flatten(P))
partition_flat = collect(Iterators.flatten(partition))

# P should contain unique integers from 1:length(X)
@assert norm(unique(sort(P_flat)) - collect(1:length(X))) < eps(T) * 10


# visualize the distribution of maximum deviations between the two options.
@show sort(max_ds0) # option 1. More parts. Many singleton parts.
@show sort(max_ds) # option 2. Less parts, higher max_ds in general since less singleton parts.

nothing
