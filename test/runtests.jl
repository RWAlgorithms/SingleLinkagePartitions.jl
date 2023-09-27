import Random
Random.seed!(25)

using Test
using LinearAlgebra

import SingleLinkagePartitions
SL = SingleLinkagePartitions

import Statistics

# # utility functions

# function getdistances(X::Vector{Vector{T}}, metricfunc) where T
    
#     K = zeros(T, length(X), length(X))
#     for i in eachindex(X)
#         for j in eachindex(X) 
#             K[i,j] = metricfunc(X[i], X[j])
#         end
#     end

#     return K
# end

# function checktoptriangle(K::Matrix{T}, lb::T)::Bool where T

#     status_flag = true
#     for j in axes(K,2)
#         for i in Iterators.take(axes(K,1), j-1)
#             status_flag = status_flag & (K[i,j] > lb)
#         end
#     end

#     return status_flag
# end

# # Tests

@testset "merge point by upper bounding the number of parts (i.e. resultant points)" begin

    Ns = [234; 55; 341]
    Ds = [54; 335; 520]
    N_maxes = [52; 17;]
    T = Float64
    metricfunc = (xx,yy)->norm(xx-yy)

    # merge points by a maximum.
    for N in Ns
        for D in Ds
            for N_max in N_maxes
                println("Working on (N, D, N_max): ", (N,D,N_max))
                X = collect( randn(T,D) for _ = 1:N )
    
                Y, status_flag = SL.mergepoints(
                    SL.UpperboundParts(N_max, Inf),
                    X,
                    metricfunc,
                )
                @test status_flag
                @test length(Y) < N_max
            end
        end
    end
end

# the elements of a partition is called a part. Some call it a cluster, if they call the set X the data set.
# the non-zero entries/parts in dists_X that are less than `distance_threshold` are combined into the same parts.
# therefore, the non-zero entries of dists_Y should all be greater or equal to distance_threshold.
@testset "merge points by stepping through single-linkage partitions until a threshold is met" begin

    D = 79
    N_pts = 100
    N_tests = 2000

    metricfunc = (xx,yy)->norm(xx-yy)

    for _ = 1:N_tests
        X = collect( randn(D) for _ = 1:N_pts )
        #@show X # debug.
        dists_X = SL.getdistances(X, metricfunc)
        distance_threshold = rand()*maximum(dists_X)
        #@show distance_threshold # debug.

        Y, status_flag = SL.mergepoints(X, metricfunc; tol = distance_threshold)
        dists_Y = SL.getdistances(Y, metricfunc)
        
        status_flag2 = SL.checktoptriangle(dists_Y, distance_threshold)
        @test status_flag2 == status_flag
    end
end

@testset "devectorized computemean()" begin

    D = 79
    N_pts = 100
    N_tests = 2000

    tol = 1e-12

    for _ = 1:N_tests
        X = collect( randn(D) for _ = 1:N_pts )
        
        @test norm(Statistics.mean(X) - SL.computemean(X)) < tol
    end
end
