
## how we determine which level (i.e., which partition) given a tolerance.
abstract type LevelOption end

"""
```
struct UseCumulativeSLDistance{T} <: LevelOption
    atol::T
end
```
See `UseSLDistance` for a description on *merging distances*.

`UseCumulativeSLDistance` specifies the threshold on the cumulative merging distance from level `0` to `l`. For use with `picklevel(`.
"""
struct UseCumulativeSLDistance{T} <: LevelOption
    atol::T
end

"""
```
struct UseSLDistance{T} <: LevelOption
    atol::T
end
```
Terminology:
The partitions of two consecutive levels in any partition tree `pt` is such that the larger level corresponds to a larger partition that contains the partition that corresponds to the smaller level.
The part/cluster that merged as we traverse from level ``l`` to level ``l+1`` is the *merging distance* of the two levels.

Description:
`UseSLDistance` specifies the threshold merging distance of level `l` and `l+1`. For use with `picklevel(`.
"""
struct UseSLDistance{T} <: LevelOption
    atol::T
end

"""
```
picklevel(C::UseSLDistance, pt::PartitionTree{T}, args...)::Int where T <: AbstractFloat
```
See `UseSLDistance` for a description of *merging distance*.
Returns the partition for first level `l-1` such that the merging distance between levels ``l`` and ``l-1`` is larger than the `C.atol` struct field. The search starts from ``l`` = 1, ``l-1`` = 0 (the leaf level). 
"""
function picklevel(C::UseSLDistance, pt::PartitionTree{T}, args...)::Int where T <: AbstractFloat
    
    w = pt.w
    
    ind = findfirst(xx->xx>C.atol, w)
    level = length(w)
    if !isnothing(ind)
        level = ind
    end

    level -= 1
    level = clamp(level, 0, getNedges(pt))

    return level
end

"""
```
picklevel(C::UseCumulativeSLDistance, pt::PartitionTree{T}, args...)::Int where T <: AbstractFloat
```
See `UseSLDistance` for a description of *merging distance*.
Returns the partition for level `l-1` such that the cumulative sum of merging distances from the level ``0`` (the leaf level) to level ``l``is larger than the `C.atol` struct field. The search starts from ``l`` = 1, ``l-1`` = 0. 
"""
function picklevel(C::UseCumulativeSLDistance, pt::PartitionTree{T}, args...)::Int where T <: AbstractFloat
    
    w = pt.w
    
    ind = findfirst(xx->xx>C.atol, cumsum(w))
    level = length(w)
    if !isnothing(ind)
        level = ind
    end

    level -= 1
    level = clamp(level, 0, getNedges(pt))

    return level
end

## how we pick a representative point, called the center, for each part in a given partition.

abstract type PartRepOption end

struct UseMean <: PartRepOption end # use the mean of the part elements as the representative.
struct UseProximitytoMean <: PartRepOption end # pick the element in the part that is closest to the mean of all elements in the part.

abstract type ScoreTrait end
struct UseMinimum <: ScoreTrait end # pic the part representation based on a supplied score. Pick the minimum score among the pts in the part.
struct UseMaximum <: ScoreTrait end # Similar as above, but pick the maximum score.

struct UseScore{T <: AbstractFloat, ST <: ScoreTrait} <: PartRepOption
    trait::ST
    score::Vector{T}
end


# The output is a copy of X. This is because Statistics.mean creates a copy, even if it is iterated over a single element.
# mean and variance for the points inside a part, for each part in the partition.
function getcenters(
    ::UseMean,
    X::Vector{XT},
    partition::Vector{Vector{Int}},
    )::Vector{XT} where XT <: Union{Number, Vector} # XT is T or Vector{T}, T <: Number.

    return collect(
        Statistics.mean( X[n] for n in partition[k] )
        for k in eachindex(partition)
    )
end

function getcenters(
    ::UseProximitytoMean,
    X::Vector{XT},
    partition::Vector{Vector{Int}},
    )::Vector{XT} where XT <: Union{Number, Vector} # XT is T or Vector{T}, T <: Number.

    Xm = getcenters(UseMean(), X, partition)

    centers = Vector{XT}(undef, length(partition))
    for k in eachindex(partition)
        
        _, part_ind = findmin( norm(X[n] - Xm[k]) for n in partition[k] )
        ind = partition[k][part_ind]
        
        centers[k] = X[ind]
    end

    return centers
end

# the output is a copy, not reference. This is to keep the behavior consist with the other trait/options.
function getcenters(
    C::UseScore,
    X::Vector{XT},
    partition::Vector{Vector{Int}},
    )::Vector{XT} where XT <: Union{Number, Vector} # XT is T or Vector{T}, T <: Number.

    #w = C.score
    @assert length(X) == length(C.score)
    
    centers = Vector{XT}(undef, length(partition))
    for k in eachindex(partition)

        part_ind = pickrepresentation(C, partition[k])
        ind = partition[k][part_ind]
        
        centers[k] = copy(X[ind])
    end

    return centers
end

function pickrepresentation(C::UseScore{T,UseMaximum}, part) where T
    _, part_ind = findmax( C.score[n] for n in part )
    
    return part_ind
end

function pickrepresentation(C::UseScore{T,UseMinimum}, part) where T
    _, part_ind = findmin( C.score[n] for n in part )
    
    return part_ind
end

## how we compute the empirical variance (division by n not n-1) for each part in the given partition.

# this version is for when X[n] is a real-valued or complex-valued vector.
function getvariances(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    partition::Vector{Vector{Int}},
    )::Vector{T} where T <: AbstractFloat

    @assert !isempty(X)
    D = length(X[begin])

    vs = zeros(T, length(partition))
    for k in eachindex(partition)

        for d = 1:D
            itr = (X[n][d] for n in partition[k])
            vs[k] += Statistics.var(itr; corrected = false)
        end
    end

    return vs
end

# each dimension gets its own variance estimate.
function getseparatevariances(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    partition::Vector{Vector{Int}},
    )::Vector{Vector{T}} where T <: AbstractFloat

    @assert !isempty(X)
    D = length(X[begin])

    vs = Vector{Vector{T}}(undef, length(partition))
    for k in eachindex(partition)

        vs[k] = zeros(T, D)
        for d = 1:D
            itr = (X[n][d] for n in partition[k])
            vs[k][d] = Statistics.var(itr; corrected = false)
        end
    end

    return vs
end

# this version is for when X[n] is scalar.
function getvariances(
    X::Union{Vector{T}, Vector{Complex{T}}},
    partition::Vector{Vector{Int}},
    )::Vector{T} where T <: AbstractFloat

    @assert !isempty(X)

    vs = zeros(T, length(partition))
    for k in eachindex(partition)

        itr = (X[n] for n in partition[k])
        vs[k] += Statistics.var(itr; corrected = false)
    end

    return vs
end


# output is a copy.
# the tolerances are based on Euclidean norm, not metric.
function getreduceptspartition(
    level_trait::LevelOption,
    metric::MetricType,
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    )::Vector{Vector{Int}} where T <: AbstractFloat

    pt = computesl(metric, X)
    level = picklevel(level_trait, pt, X)
    partition = getpartition(pt, level)

    return partition
end

### frontends


function reducepts(
    level_trait::LevelOption,
    center_trait::PartRepOption,
    metric::MetricType,
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    ) where T <: AbstractFloat
    
    if length(X) == 1
        return deepcopy(X), zeros(T, 1), collect( ones(Int,1) for _ = 1:1 )
    end

    partition = getreduceptspartition(level_trait, metric, X)

    centers = getcenters(center_trait, X, partition)
    variances = getvariances(X, partition)

    return centers, variances, partition
end

# the version that also also reduces y.
# usage: reducing training set of scalar-valued (X, y) regression data.
function reducepts(
    level_trait::LevelOption,
    center_trait::PartRepOption,
    metric::MetricType,
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    y::Vector{YT},
    ) where {T <: AbstractFloat, YT <: Number}
    
    if length(X) == 1
        return deepcopy(X), zeros(T, 1), copy(y), zeros(T, 1), collect( ones(Int,1) for _ = 1:1 )
    end

    partition = getreduceptspartition(level_trait, metric, X)

    Xc = getcenters(center_trait, X, partition)
    Xv = getvariances(X, partition)

    yc = getcenters(center_trait, y, partition)
    yv = getvariances(y, partition)

    return Xc, Xv, yc, yv, partition
end

# usage: reducing training set of vector-valued (X, y_set) regression data.
function reducepts(
    level_trait::LevelOption,
    center_trait::PartRepOption,
    metric::MetricType,
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    y_set::Vector{Vector{YT}},
    ) where {T <: AbstractFloat, YT <: Number}
    
    if length(X) == 1
        return deepcopy(X), zeros(T, 1), deepcopy(y_set), collect( zeros(T, 1) for _ in eachindex(y_set) ), collect( ones(Int,1) for _ = 1:1 )
    end

    partition = getreduceptspartition(level_trait, metric, X)

    Xc = getcenters(center_trait, X, partition)
    Xv = getvariances(X, partition)

    yc_set = collect( getcenters(center_trait, y, partition) for y in y_set )
    yv_set = collect( getvariances(y, partition) for y in y_set )

    return Xc, Xv, yc_set, yv_set, partition
end