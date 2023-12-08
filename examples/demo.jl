
#using LinearAlgebra

import Random
Random.seed!(25)


# import SingleLinkagePartitions
# SL = SingleLinkagePartitions


#### singlet-linkage clustering.
T = Float64

# with non-unique values.
X, ref_parts = deserialize("data/leucine700")
@assert length(collect(Iterators.flatten(ref_parts))) == length(X)

# from the README.md.
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

#X = collect( ones(T, 4) for _ = 1:11) # identical set.

#X = collect( ones(T, 4) for _ = 1:11) # identical but one.
#push!(X, zeros(4))

# t = LinRange(-0.05, 0.1, 10000)
# X = collect( [t_i;] for t_i in t)


N = length(X)

metric = SL.geteuclideanmetric()
pt = SL.computesl(metric, X)

level = 2
partition = SL.getpartition(pt, level)

distance_set = SL.getdistances(pt)
partition_set = SL.generateallpartitions(pt)

#### reduce points.

zero_tol = convert(T, eps(T)*100)
atol = convert(T, 1e-4)


score = rand(T, N)


level_trait = SL.UseSLDistance()

center_trait = SL.UseScore(SL.UseMinimum(), score) # when Xc must be a subset of X
center_trait = SL.UseProximitytoMean() # when Xc must be a subset of X
center_trait = SL.UseMean() # when Xc can be any pt from which X is constructed from.


Xc0, vs0, partition_r = SL.reducepts(level_trait, center_trait, metric, X, atol)

y = randn(Complex{Float32}, N)
Xc1, vs_X1, yc, vs_y, partition_r = SL.reducepts(
    level_trait,
    center_trait,
    metric,
    X,
    y,
    atol,
)


y_set = collect( randn(Complex{Float32}, N) for _ = 1:3 )
Xc, vs_X, yc_set, vs_y_set, partition_r = SL.reducepts(
    level_trait,
    center_trait,
    metric,
    X,
    y_set,
    atol,
)


@assert norm(Xc0 - Xc1) < eps(T)*100
@assert norm(Xc - Xc1) < eps(T)*100

@assert norm(vs0 - vs_X1) < eps(T)*100
@assert norm(vs_X - vs_X1) < eps(T)*100

########### permutation test.


level_trait = SL.UseSLDistance()
level_trait = SL.UseCumulativeSLDistance()

X0 = collect( randn(T, 3) for _ = 1:9 )
X0 = collect( randn(T,3) .* atol/4 for _ = 1:1000 )
Xc0, vs0, partition0 = SL.reducepts(level_trait, center_trait, metric, X0, atol)

X1 = deepcopy(X0)
X1[1], X[end] = X[end], X[1]
Xc1, vs1, partition1 = SL.reducepts(level_trait, center_trait, metric, X0, atol)

@show norm(sort(Xc0) - sort(Xc1))


nothing
