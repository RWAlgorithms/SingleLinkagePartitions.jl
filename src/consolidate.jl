
# how to consolidate the points in a part to one point.
abstract type ConsolidationOption end

struct UseMinimum <: ConsolidationOption end
struct UseMaximum <: ConsolidationOption end

# returns k or less distinct points, based on the trait and the values in y.
function reducepoints(
    trait::ConsolidationOption,
    k::Int,
    X::Vector{Vector{T}},
    y::Vector{T},
    metric::MetricType;
    early_stop_distance = convert(T, Inf),
    zero_tol::T = eps(T)*100,
    ) where T

    @assert length(X) == length(y)
    @assert k > 0
    k = min(k, length(X))

    distance_set, partition_set = runsinglelinkage(
        X,
        metric;
        early_stop_distance = early_stop_distance
    )
    
    # allocate for scope purposes.
    partition = Vector{Vector{Int}}(undef, 0)

    if isapprox(distance_set[end], zero(T); atol = zero_tol)
        # all the same points in X. take the y that satisfies the trait.
        partition = partition_set[end]
    else

        # the partitions including and after partition_set[ind_1st_unique] have unique points from X in each part.
        ind_1st_unique = length(partition_set)
        result = findfirst(xx->xx>zero_tol, distance_set)
        if !isnothing(result)
            ind_1st_unique = result
        end

        ind = length(partition_set) - k + 1

        # choose the partition level such that Xr contain distinct points.
        if ind < ind_1st_unique
            ind = ind_1st_unique
        end

        #@show result, ind, partition_set
        partition = partition_set[ind]
        # if length(partition) != k
        #     partition, _ = findmin(xx->abs(length(xx)-k), partition_set)
        # end
    end

    Xr, yr = consolidate(trait, X, y, partition)
    
    return Xr, yr, partition
end

function consolidate(::UseMinimum, X::Vector{T2}, y::Vector{T}, partition::Vector{Vector{Int}}) where {T,T2}
    
    out_X = Vector{T2}(undef, length(partition))
    out_y = Vector{T}(undef, length(partition))
    for i in eachindex(partition)

        _, ind = findmin( y[m] for m in partition[i] )
        out_X[i] = X[partition[i][ind]]
        out_y[i] = y[partition[i][ind]]
    end

    return out_X, out_y
end

function consolidate(::UseMaximum, X::Vector{T2}, y::Vector{T}, partition::Vector{Vector{Int}}) where {T,T2}
    
    out_X = Vector{T2}(undef, length(partition))
    out_y = Vector{T}(undef, length(partition))
    for i in eachindex(partition)

        _, ind = findmax( y[m] for m in partition[i] )
        out_X[i] = X[partition[i][ind]]
        out_y[i] = y[partition[i][ind]]
    end

    return out_X, out_y
end


#### the new mergepoints.

function reducepoints(
    k::Int,
    X::Vector{Vector{T}},
    metric::MetricType;
    early_stop_distance = convert(T, Inf),
    zero_tol::T = eps(T)*100,
    ) where T

    @assert k > 0
    k = min(k, length(X))

    distance_set, partition_set = runsinglelinkage(
        X,
        metric;
        early_stop_distance = early_stop_distance
    )
    
    # allocate for scope purposes.
    partition = Vector{Vector{Int}}(undef, 0)

    if isapprox(distance_set[end], zero(T); atol = zero_tol)
        # all the same points in X.
        partition = partition_set[end]
    else

        # the partitions including and after partition_set[ind_1st_unique] have unique points from X in each part.
        ind_1st_unique = length(partition_set)
        result = findfirst(xx->xx>zero_tol, distance_set)
        if !isnothing(result)
            ind_1st_unique = result
        end

        ind = length(partition_set) - k + 1

        # choose the partition level such that Xr contain distinct points.
        if ind < ind_1st_unique
            ind = ind_1st_unique
        end

        #@show result, ind, partition_set
        partition = partition_set[ind]
        # if length(partition) != k
        #     partition, _ = findmin(xx->abs(length(xx)-k), partition_set)
        # end
    end

    Xr = averagepoints(X, partition)
    
    return Xr, partition
end

function averagepoints(X::Vector{T}, partition::Vector{Vector{Int}}) where T
    
    Y = Vector{T}(undef, length(partition))
    for i in eachindex(partition)

        Y[i] = Statistics.mean( X[m] for m in partition[i] )
    end

    return Y
end