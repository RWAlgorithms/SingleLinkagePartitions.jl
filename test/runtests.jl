using Test, LinearAlgebra, Random, Statistics
Random.seed!(25)

import SingleLinkagePartitions as SL

function is_same_partition_tree(pt1::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}}, pt2::AbstractVector{<:AbstractVector{<:AbstractVector{<:Integer}}})

    flag = true
    for (u, U) in Iterators.zip(pt1, pt2)

        flag = flag && all(x in U for x in u)
    end

    return flag
end

@testset "SLINK test 1" begin

    T = Float64
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

    partition_set_oracle = [[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10]], [[1], [2], [3], [4], [5], [6], [7, 10], [8], [9]], [[1, 9], [2], [3], [4], [5], [6], [7, 10], [8]], [[1, 9], [2], [3], [4], [5], [6, 7, 10], [8]], [[1, 9], [2], [3], [4, 6, 7, 10], [5], [8]], [[1, 2, 9], [3], [4, 6, 7, 10], [5], [8]], [[1, 2, 3, 9], [4, 6, 7, 10], [5], [8]], [[1, 2, 3, 8, 9], [4, 6, 7, 10], [5]], [[1, 2, 3, 8, 9], [4, 5, 6, 7, 10]], [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]]
    distance_set_oracle = [0.6502879922496316, 0.7069552456099949, 0.7140482053065044, 0.7164233686190572, 0.7855224671411318, 0.9225136028232729, 1.2606479889562476, 1.4987127483818177, 1.5900454258666015]

    # # SLINK
    N = length(X)
    X_mat = reshape(collect(Iterators.flatten(X)), length(X[begin]), length(X))

    # allocate
    s = SL.SLINKState(T, N)
    dist_callable = SL.EuclideanDistance()

    # Compute the essential data required to create a partition tree. Non-allocating.
    SL.slink!(s, dist_callable, X_mat)

    @test norm(SL.construct_distances(s) - distance_set_oracle) < eps(T) * 100

    # construct the partition tree. Allocates.
    partition_set = SL.construct_partition_tree(s)
    @test is_same_partition_tree(partition_set, partition_set_oracle)

    # ## Select the partition from the tree that corresponds to having had 6 merges.
    # num_fuses can take values from 0 (singleton parts) to length(X) - 1 (one single part).
    num_fuses = 6
    selected_partition = SL.construct_partition(s, num_fuses)

    oracle = Memory{Int}[[5], [8], [1, 2, 3, 9], [4, 6, 7, 10]]
    @test all(x in oracle for x in selected_partition)
end
