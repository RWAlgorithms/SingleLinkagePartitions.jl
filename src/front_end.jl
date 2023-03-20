
##### applications.

function computemean(x::Vector{T})::T where T
    if isempty(x)
        return zero(T)
    end

    return sum(x)/length(x)
end

"""
function mergepoints(
    X::Vector{Vector{T}},
    metricfunc::Function;
    tol = 1e-6,
    )::Tuple{Vector{Vector{T}},Bool} where T

Inputs:
- `X` is the set of points to be merged. This algorithm assumes the points are unique.

- `metricfunc` is the metric function used to compute the distance between two points.

- `tol` is the minimum pair-wise distance we desire for the output set of points, `Y`.


Outputs:
- `Y` is the set of merged points.

- `status_flag` indicates whether the pair-wise distances in `Y` are all smaller than `tol`.
"""
function mergepoints(
    X::Vector{Vector{T}},
    metricfunc::Function;
    tol = 1e-6,
    )::Tuple{Vector{Vector{T}},Bool} where T

    @assert tol > zero(T)

    h_set, partition_set = runsinglelinkage(X, metricfunc; early_stop_distance = tol)

    inds = sortperm(h_set, rev = true)

    partitioned_set_sorted = partition_set[inds]
    h_set_sorted = h_set[inds]

    ind = findfirst(xx->xx<tol, h_set_sorted)
    if typeof(ind) <: Nothing
        println("Error with mergepoints(). Returning a copy of the input.")
        return copy(X)
    end

    # Since we're replacing each part with its centroid, We might violate the tolerance.
    Y = Vector{Vector{T}}(undef, 0) # allocate for scope reasons.
    tol_satisfied = false

    while !tol_satisfied && ind > 0
        partition = partitioned_set_sorted[ind] # select partition.
        Y = fuseparts(partition, X)

        Y_dists = getdistances(Y, metricfunc)
        tol_satisfied = checktoptriangle(Y_dists, tol)

        ind -= 1
    end
    
    return Y, tol_satisfied
end

function fuseparts(
    partition::Vector{Vector{Int}},
    X::Vector{Vector{T}},
    ) where T

    N_parts = length(partition)
    Y = Vector{Vector{T}}(undef, N_parts)
    #max_distance_from_tol = Vector{T}(undef, N_parts)
    for k in eachindex(Y)
        
        S = X[partition[k]]
        #out[k] = Statistics.mean(S)
        mean_k = computemean(S) # rid of dependency on Statistics.

        #max_distance_from_tol[k] = maximum( metricfunc(S[i], mean_k) for i in eachindex(S) )
        Y[k] = mean_k
    end

    return Y
end

########## routines for testing.

function getdistances(X::Vector{Vector{T}}, metricfunc) where T
    
    K = zeros(T, length(X), length(X))
    for i in eachindex(X)
        for j in eachindex(X) 
            K[i,j] = metricfunc(X[i], X[j])
        end
    end

    return K
end

function checktoptriangle(K::Matrix{T}, lb::T)::Bool where T

    status_flag = true
    for j in axes(K,2)
        for i in Iterators.take(axes(K,1), j-1)
            status_flag = status_flag && (K[i,j] > lb)
            #@show (i,j)
        end
    end

    return status_flag
end