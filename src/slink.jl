# SPDX-License-Identifier: AGPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

struct SLINKState{T <: AbstractFloat}
    M::Memory{T} # buffer
    L::Memory{T} # height values
    P::Memory{Int} # last object index in some part that just fused.

    @doc """
        SLINKState(::Type{T}, N::Integer) where {T <: AbstractFloat}

    Instantiate a new `SLINKState{T}` buffer for `N` observations. The internals of `SLINKState` are considered private API.
    """
    function SLINKState(::Type{T}, N::Integer) where {T <: AbstractFloat}
        return new{T}(Memory{T}(undef, N), Memory{T}(undef, N), Memory{Int}(undef, N))
    end
end

function get_dissimilarities(s::SLINKState)
    return s.M
end

"""
    get_cardinality(s::SLINKState)

Returns the number of entries for the input set. Return type is `Int`.
"""
function get_cardinality(s::SLINKState)
    return length(s.P)
end

"""
    slink!(s::SLINKState{T}, eval_dist, X::AbstractMatrix) where {T <: AbstractFloat}

`X` is the point cloud matrix, and must be 1-indexing with stride 1. `X[d,n]` should be the `d`-th coordinate dimension, `n`-th data entry.
Returns `nothing`.

`eval_dist` should be a callable (struct or function) that takes inputs `(x::AbstractVector{T}, y::AbstractVector{T})` and outputs a number of type `T`. An example is `EuclideanDistance`.
"""
function slink!(s::SLINKState{T}, eval_dist, X::AbstractMatrix) where {T <: AbstractFloat}

    M, L, P = s.M, s.L, s.P

    P[begin] = 1
    L[begin] = one(T)

    for t in Iterators.take(eachindex(L, P), length(L) - 1)

        # prepare
        P[t + 1] = t + 1
        L[t + 1] = T(Inf)

        v = view(X, :, t + 1)
        for i in Iterators.take(eachindex(M, axes(X, 2)), t)
            M[i] = eval_dist(view(X, :, i), v)
        end

        #
        for i in Iterators.take(eachindex(M, L, P), t)
            if L[i] >= M[i]
                M[P[i]] = min(M[P[i]], L[i])
                L[i] = M[i]
                P[i] = t + 1
            else
                M[P[i]] = min(M[P[i]], M[i])
            end
        end

        #
        for i in Iterators.take(eachindex(L, P), t)
            if L[i] >= L[P[i]]
                P[i] = t + 1
            end
        end
    end

    return nothing
end

"""
    slink!(s::SLINKState{T}, eval_dist, X::AbstractVector{<:AbstractVector}) where {T <: AbstractFloat}
"""
function slink!(s::SLINKState{T}, eval_dist, X::AbstractVector{<:AbstractVector}) where {T <: AbstractFloat}

    M, L, P = s.M, s.L, s.P

    P[begin] = 1
    L[begin] = one(T)

    for t in Iterators.take(eachindex(L, P), length(L) - 1)

        # prepare
        P[t + 1] = t + 1
        L[t + 1] = T(Inf)

        v = X[t + 1]
        for i in Iterators.take(eachindex(M, X), t)
            M[i] = eval_dist(X[i], v)
        end

        #
        for i in Iterators.take(eachindex(M, L, P), t)
            if L[i] >= M[i]
                M[P[i]] = min(M[P[i]], L[i])
                L[i] = M[i]
                P[i] = t + 1
            else
                M[P[i]] = min(M[P[i]], M[i])
            end
        end

        #
        for i in Iterators.take(eachindex(L, P), t)
            if L[i] >= L[P[i]]
                P[i] = t + 1
            end
        end
    end

    return nothing
end

"""
    struct EuclideanDistance end

Callable struct:
```
(::EuclideanDistance)(x::AbstractVector, y::AbstractVector)
```
Returns the Euclidean distance between `x` and `y`, without allocating.
"""
struct EuclideanDistance end
function (::EuclideanDistance)(x::AbstractVector, y::AbstractVector)
    return sqrt(sum((x[i] - y[i])^2 for i in eachindex(x, y)))
end

#### output processing

"""
    construct_distances(s::SLINKState)

Returns the set of minimum distanes between consecutive nested partitions in the partition tree associated with `s`. Run `slink!` on `s` first.

This function allocates. The returned array should have `get_cardinality(s)-1` entries. It is of type `Memory{T}` if `s::SLINKState{T}`.
"""
function construct_distances(s::SLINKState)
    return sort(s.L)[1:(length(s.L) - 1)]
end

"""
    construct_partition_tree(s::SLINKState)

This function allocates. Run `slink!` on `s` first. Return type is `Memory{Memory{Memory{Int}}}`.
"""
function construct_partition_tree(s::SLINKState)

    out = Memory{Memory{Memory{Int}}}(undef, get_cardinality(s))
    for i in eachindex(out)
        k = i - 1
        out[i] = construct_partition(s, k)
    end
    return out
end

"""
    construct_partition(s::SLINKState, num_fuses::Integer)

This function allocates. Run `slink!` on `s` first. Return type is `Memory{Memory{Int}}`.

`num_fuses` identifies the partition tree level. It can take on values between 0 to `get_cardinality(s)`.
"""
function construct_partition(s::SLINKState, num_fuses::Integer)

    N = get_cardinality(s)
    0 <= num_fuses < N || error("num_fuses is out of bounds.")

    # trivial cases.
    partition = Memory{Memory{Int}}(undef, N)
    for i in eachindex(partition)
        partition[i] = Memory{Int}(undef, 1)
        partition[i][begin] = i
    end

    if num_fuses == 0
        return partition
    end

    if num_fuses >= N - 1
        partition[1] = collect(1:N)
        return _post_process_partition(view(partition, 1:1))
    end

    # general cases.
    fuse_inds = sortperm(s.L) # allocates
    processed_num_fuses = 0
    for i in fuse_inds

        if processed_num_fuses >= num_fuses
            return _post_process_partition(partition)
        end
        processed_num_fuses += 1

        # move object i to the part that has object P[i]
        src_part_ind = _find_part(partition, i)
        dest_part_ind = _find_part(partition, s.P[i])

        # merge.
        tmp = vcat(partition[src_part_ind], partition[dest_part_ind])
        sort!(tmp)
        partition[dest_part_ind] = tmp
        partition[src_part_ind] = Memory{Int}(undef, 0)
    end

    return _post_process_partition(partition)
end

# Each part must be sorted in ascending order (the contents are object indices).
function _find_part(partition::AbstractVector{<:AbstractVector{<:Integer}}, k::Integer)

    for i in eachindex(partition)
        if !isempty(partition[i])
            if partition[i][end] == k
                return i
            end
        end
    end
    return 0 # error case.
end

function _post_process_partition(partition)

    num_parts = length(partition) - sum(isempty(p) for p in partition)
    out = Memory{Memory{Int}}(undef, num_parts)

    k = 1
    for i in eachindex(partition)
        if !isempty(partition[i])
            out[k] = Memory{Int}(partition[i])
            k += 1
        end
    end

    return out
end
