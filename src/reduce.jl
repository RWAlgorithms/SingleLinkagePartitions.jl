
## how we determine which level (i.e., which partition) given a tolerance.
abstract type LevelOption end

# the output should have a lower minimum pair-wise distance than atol, but might have give different results when the input set X is permutated.
struct UseCumulativeSLDistance{T} <: LevelOption
    atol::T
end

# better consistency of algorithm, but the output might have a minimum pair-wise distance larger than the given atol.
struct UseSLDistance{T} <: LevelOption
    atol::T
end

function getlevel(C::UseSLDistance, pt::PartitionTree{T}, args...)::Int where T <: AbstractFloat
    
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

function getlevel(C::UseCumulativeSLDistance, pt::PartitionTree{T}, args...)::Int where T <: AbstractFloat
    
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
    level = getlevel(level_trait, pt, X)
    partition = getpartition(pt, level)

    return partition
end

### frontends


function reducepts(
    level_trait::LevelOption,
    center_trait::PartRepOption,
    metric::MetricType,
    X::Vector{Vector{XT}},
    ) where XT <: Number
    
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
    X::Vector{Vector{XT}},
    y::Vector{YT},
    ) where {XT <: Number, YT <: Number}
    
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
    X::Vector{Vector{XT}},
    y_set::Vector{Vector{YT}},
    ) where {XT <: Number, YT <: Number}
    
    partition = getreduceptspartition(level_trait, metric, X)

    Xc = getcenters(center_trait, X, partition)
    Xv = getvariances(X, partition)

    yc_set = collect( getcenters(center_trait, y, partition) for y in y_set )
    yv_set = collect( getvariances(y, partition) for y in y_set )

    return Xc, Xv, yc_set, yv_set, partition
end