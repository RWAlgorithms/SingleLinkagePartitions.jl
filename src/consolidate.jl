
# how to consolidate the points in a part to one point.
abstract type ConsolidationOption end

struct UseMinimum <: ConsolidationOption end
struct UseMaximum <: ConsolidationOption end

function reducepoints(
    trait::ConsolidationOption,
    k::Int,
    X::Vector{Vector{T}},
    y::Vector{T},
    metric::MetricType;
    early_stop_distance = convert(T, Inf),
    ) where T

    _, partition_set = runsinglelinkage(
        X,
        metric;
        early_stop_distance = early_stop_distance
    )
    
    partition = partition_set[end-k+1]
    if length(partition) != k
        partition, _ = findmin(xx->abs(length(xx)-k), partition_set)
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
