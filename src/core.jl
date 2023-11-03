
# order of points in X matters.
# based on https://sites.cs.ucsb.edu/~veronika/MAE/summary_SLINK_Sibson72.pdf



# # debug data structure.
#= abstract type StoreDiagnosticsTrait end

struct StoreDiagnostics <: StoreDiagnosticsTrait end
struct DiscardDiagnostics <: StoreDiagnosticsTrait end

struct LevelDiagnostics
    i::Int
    j::Int
    
    k1::Int
    k2::Int

    map_part::Vector{Int}
    #partition::Vector{Vector{Int}}
end

function storelevel!(A::Vector{LevelDiagnostics}, ::StoreDiagnostics, i, j, k1, k2, R)
    push!(
        A,
        LevelDiagnostics(i, j, k1, k2, copy(R.map_part))
    )
    return nothing
end

function storelevel!(A::Vector{LevelDiagnostics}, ::DiscardDiagnostics, args...)
    return nothing
end =#


#### main data structures.

struct DistanceContainer{T}
    buffer::Matrix{T}

    # indices for distance parts for parts in the current partition.
    map_part::Vector{Int} # if i := map_part[k], then part k corresponds to row or column i in the distance buffer.
    partition::Vector{Vector{Int}} # current parition of point indices.

    #rows::Vector{Int} # i := rows[k] is the row index for the k-th part in the current partition.
    #cols::Vector{Int} # j := cols[k] is the col index for the k-th part in the current partition.
end

function DistanceContainer(X::Vector{Vector{T}}, metric::MetricType)::DistanceContainer{T} where T

    N = length(X)
    A = Matrix{T}(undef, N,N)
    for i2 in axes(A, 2)
        for i1 in axes(A, 1)
            A[i1,i2] = evalmetric(metric, X[i1], X[i2])
        end
    end

    #R = DistanceContainer(A, collect(1:N), collect(1:N))
    R = DistanceContainer(
        A,
        collect(1:N),
        collect([i;] for i in eachindex(X)),
    )
    postprocessbuffer!(R)

    return R
end

function postprocessbuffer!(B::DistanceContainer{T}) where T
    
    A = B.buffer

    # exclude the distance entries of a point with itself.
    for i in axes(A,1)
        A[i,i] = convert(T, Inf)
    end
    
    return nothing
end

function isvalidcontainer(A::DistanceContainer)::Bool
    return !isempty(A.map_part)
end

function getNparts(A::DistanceContainer)::Int
    return length(A.partition)
end

# (i,j) are buffer indices. k1, k2 are part indices.
function findmindistance(A::DistanceContainer{T}) where T
    
    @assert !isempty(A.buffer)
    #@assert !isempty(A.rows)
    #@assert !isempty(A.cols)
    @assert !isempty(A.map_part)

    min_i = 0
    min_j = 0
    min_k1 = 0
    min_k2 = 0
    min_Aij = convert(T, Inf)
    
    #for k2 in eachindex(A.cols)
    for k2 in eachindex(A.map_part)
        j = A.map_part[k2]

        #for k1 in eachindex(A.rows)
        for k1 in eachindex(A.map_part)
            i = A.map_part[k1]

            if A.buffer[i,j] < min_Aij
                
                min_i, min_j, min_k1, min_k2 = i, j, k1, k2
                
                min_Aij = A.buffer[i,j]
            end
        end
    end

    return min_i, min_j, min_k1, min_k2, min_Aij
end

# merge two clusters.
# (i0,j0) are buffer indices. k1, k2 are part indices to be merged.
function mergeparts!(
    A::DistanceContainer{T},
    i0::Integer,
    j0::Integer,
    k1::Integer,
    k2::Integer,
    ) where T

    # update part1 with min distance,
    #for k in eachindex(A.cols)
    for k in eachindex(A.map_part)
        j = A.map_part[k]

        A.buffer[i0, j] = min(
            A.buffer[i0, j],
            A.buffer[j0, j],
        )
    end

    # since we're working without symmetry in findminstance(), we need to update the redundeant lower triangular half too.
    #for k in eachindex(A.rows)
    for k in eachindex(A.map_part)
        i = A.map_part[k]

        A.buffer[i, i0] = min(
            A.buffer[i, i0],
            A.buffer[i, j0],
        )
    end

    # make it such that the diagonals won't be returned in findmindistance()
    A.buffer[j0,j0] = convert(T, Inf)
    A.buffer[i0,i0] = convert(T, Inf)
    #postprocessbuffer!(dist)

    # merge k1 and k2 to the same part.
    A.partition[k1] = [ A.partition[k1]; A.partition[k2]; ]
    
    # remove part2 from search list (aka. active indices)
    #deleteat!(A.rows, k2)
    #deleteat!(A.cols, k2)
    deleteat!(A.map_part, k2)
    deleteat!(A.partition, k2)

    return nothing
end

# for compatibility with legacy code.
function runsinglelinkage(
    X::Vector{Vector{T}},
    metricfunc::Function;
    early_stop_distance = convert(T, Inf), # default to solve for all possible levels.
    )::Tuple{Vector{T},Vector{Vector{Vector{Int}}}} where T

    return runsinglelinkage(X, GeneralMetric(metricfunc); early_stop_distance = early_stop_distance)
end

"""
function runsinglelinkage(
    X::Vector{Vector{T}},
    metric::MetricType;
    early_stop_distance = Inf, # default to solve for all possible levels.
    )::Tuple{Vector{T},Vector{Vector{Vector{Int}}}} where T

Creates a sequence of partitions of the points in `X`. The previous partition is nested in the current partition, and only one part is merged in each iteration of this partition building.
The two parts that merge at an iteration, or called level, is the two parts that have the smallest pair-wise distance in the partition for that level.

We use the terminology `part` for an element of a partition.

Inputs:
- `X` is the set of points to be merged. This algorithm assumes the points are unique.

- `metric` is the metric function used to compute the distance between two points.

- `early_stop_distance`: `runsinglelinkage()` terminates early when the minimum pair-wise distance of the current partition is less than `early_stop_distance`.
Setting `early_stop_distance = Inf` tells `runsinglelinkage()` to keep building partitions until the current partition has only one part, which contains all the points in `X`.


Outputs:
- `h_set::Vector{T}`: its i-th entry contains the distances of the parts that were merged at the i-th level. The distance for the first level is always zero.

- `partition_set` its i-th entry contains the partition (expressed in arrays of indices for `X`) at level i. The partition for the first level is always the point where every point in `X` is a part.

"""
function runsinglelinkage(
    X::Vector{Vector{T}},
    metric::MetricType;
    early_stop_distance = convert(T, Inf), # default to solve for all possible levels.
    )::Tuple{Vector{T},Vector{Vector{Vector{Int}}}} where T

    # set up dendrogram.
    h_set = Vector{T}(undef,0) # fusion distance at each level.
    push!(h_set, zero(T))

    partition_set = Vector{Vector{Vector{Int}}}(undef,0)
    push!(partition_set, collect([i;] for i in eachindex(X)))
    
    #diagnostics = Vector{LevelDiagnostics}(undef, 0)

    # intermediates.
    R = DistanceContainer(X, metric)

    # subsequent levels.
    #new_part_ind = length(X) # we only merge two parts at each level, so only one new partition index per level.

    # Each iteration creates a level in the dendrogram.
    # stop if:
    # - No more parts, for some reason. This is an unexpected error case, because we should'ave had the stop at 1 partition case first.
    # - The partition is a singleton partition; i.e. getNparts(R) == 1.
    # - The min-dist used to form the current partition is less than early_stop_distance.
    while isvalidcontainer(R) && getNparts(R) > 1 && h_set[end] < early_stop_distance
        
        i, j, k1, k2, min_dist = findmindistance(R)
        
        # make this into diagnostic so we can visualize.
        #@show i, j, k1, k2, min_dist, R.map_part, getNparts(R) # debug.
        #storelevel!(diagnostics, store_trait, i, j, k1, k2, R)

        mergeparts!(R, i, j, k1, k2)

        push!(partition_set, deepcopy(R.partition)) # TODO: think of an efficient method to do this without deep copying. maybe keep track of the difference instead of concatenating k1, k2, in merparts!().
        push!(h_set, min_dist)
    end

    #return h_set, partition_set, diagnostics
    return h_set, partition_set
end