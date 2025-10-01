# SPDX-License-Identifier: AGPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

# If a representative is a member of the input point set, then user can choose to make a copy or return a reference to the member.
abstract type CopyTrait end
struct NoCopy <: CopyTrait end
function package_pt(::NoCopy, x::Union{AbstractVector, Number})
    return x
end

struct MakeCopy <: CopyTrait end # make a copy of the return point set.
function package_pt(::MakeCopy, x::Union{AbstractVector, Number})
    return copy(x)
end

## how we pick a representative point, called the center, for each part in a given partition.

"""
    abstract type PartRepOption end

A trait that is a subtype of `PartRepOption` dictates how `reduce_pts` select the reduced point output from all the points of a given part that is to be reduced. Use `subtypes` to see the implemented traits.
"""
abstract type PartRepOption end

#
"""
    struct UseMean <: PartRepOption end

This is a trait that specifies `reduce_pts` should use the mean of the part elements as the reduced point. This means the output is in general not an element of the part to be reduced.
"""
struct UseMean <: PartRepOption end

#
"""
    struct UseProximitytoMean <: PartRepOption end

This is a trait that specifies `reduce_pts` should use the element closest to the mean of all the elements in the part to the reduced.
"""
struct UseProximitytoMean{CT <: CopyTrait} <: PartRepOption
    copy_trait::CT
end

function package_pt(C::UseProximitytoMean, x)
    return package_pt(C.copy_trait, x)
end

#
"""
    abstract type ScoreTrait end

A trait that is a subtype of `ScoreTrait` dictates how `reduce_pts` uses an extern score array, with eleemnts corresponding to the input point cloud that generated the partition tree, to select the reduced point output from all the points of a given part that is to be reduced. Use `subtypes` to see the implemented traits.
"""
abstract type ScoreTrait end

"""
    struct UseMinimum <: ScoreTrait end

This is a trait that specifies `reduce_pts` should use the element that corresponds to the minimum score in the score array.
"""
struct UseMinimum <: ScoreTrait end

"""
    struct UseMaximum <: ScoreTrait end

This is a trait that specifies `reduce_pts` should use the element that corresponds to the maximum score in the score array.
"""
struct UseMaximum <: ScoreTrait end # Similar as above, but pick the maximum score.

struct UseScore{T <: AbstractFloat, ST <: ScoreTrait, CT <: CopyTrait} <: PartRepOption

    score::Memory{T}
    score_trait::ST
    copy_trait::CT

    @doc """
        UseScore(
            copy_trait::CT,
            score_trait::ST,
            score::AbstractVector{T},
        ) where {T <: AbstractFloat, ST <: ScoreTrait, CT <: CopyTrait}

    Creates a variable of type `UseScore{T, ST, CT} <: PartRepOption`. The internals of `UseScore` are considered private API.
    """
    function UseScore(
            copy_trait::CT,
            score_trait::ST,
            score::AbstractVector{T},
        ) where {T <: AbstractFloat, ST <: ScoreTrait, CT <: CopyTrait}
        return new{T, ST, CT}(Memory{T}(score), score_trait, copy_trait)
    end
end

function package_pt(C::UseScore, x)
    return package_pt(C.copy_trait, x)
end

#### methods

# The output is a copy of X. This is because Statistics.mean creates a copy, even if it is iterated over a single element.
# mean and variance for the points inside a part, for each part in the partition.
# output is of type Memory{XT}, XT <: Union{Number, AbstractVector{<:Number}}
# No NoCopy option, because this cannot return a reference.
function compute_representatives(
        ::UseMean,
        X::AbstractVector{XT},
        partition::AbstractVector{<:AbstractVector{Int}},
    ) where {XT <: Union{Number, AbstractVector{<:Number}}}

    centers = Memory{XT}(undef, length(partition))
    for k in eachindex(partition, centers)
        centers[k] = Statistics.mean(X[n] for n in partition[k])
    end
    return centers
end

# For a given part where its memebers are to be merged into a single point (the center), set that single point as the member that is closest to the mean of the members.
# Returns a reference.
function compute_representatives(
        C::UseProximitytoMean,
        X::AbstractVector{XT},
        partition::AbstractVector{<:AbstractVector{Int}},
    ) where {XT <: Union{Number, AbstractVector{<:Number}}}

    Xm = compute_representatives(UseMean(), X, partition)

    centers = Memory{XT}(undef, length(partition))
    for k in eachindex(partition)

        _, part_ind = findmin(norm(X[n] - Xm[k]) for n in partition[k])
        ind = partition[k][part_ind]

        centers[k] = package_pt(C.copy_trait, X[ind])
    end
    return centers
end

# the output is a copy, not reference. This is to keep the behavior consist with the other trait/options.
# For a given part where its memebers are to be merged into a single point (the center), set that single point as the member that is closest to the mean of the members.
function compute_representatives(
        C::UseScore,
        X::AbstractVector{XT},
        partition::AbstractVector{<:AbstractVector{Int}},
    ) where {XT <: Union{Number, AbstractVector{<:Number}}}

    @assert length(X) == length(C.score)

    centers = Memory{XT}(undef, length(partition))
    for k in eachindex(partition)

        part_ind = pickrepresentation(C.score_trait, C, partition[k])
        ind = partition[k][part_ind]

        centers[k] = package_pt(C.copy_trait, X[ind])
    end
    return centers
end

function pickrepresentation(::UseMaximum, C::UseScore, part)
    _, part_ind = findmax(C.score[n] for n in part)
    return part_ind
end

function pickrepresentation(::UseMinimum, C::UseScore, part)
    _, part_ind = findmin(C.score[n] for n in part)
    return part_ind
end

## how we compute the empirical variance (division by n not n-1) for each part in the given partition.
# this version is for when X[n] is a real-valued or complex-valued vector.
function _compute_variances(
        X::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}},
        partition::AbstractVector{<:AbstractVector{Int}},
    ) where {T <: AbstractFloat}

    !isempty(X) || error("X must be non-empty.")
    D = length(X[begin])

    vs = Memory{T}(undef, length(partition))
    fill!(vs, zero(T))
    for k in eachindex(partition)

        for d in 1:D
            itr = (X[n][d] for n in partition[k])
            vs[k] += Statistics.var(itr; corrected = false)
        end
    end
    return vs
end

# this version is for when X[n] is scalar.
function _compute_variances(
        X::Union{AbstractVector{T}, AbstractVector{Complex{T}}},
        partition::AbstractVector{<:AbstractVector{Int}},
    ) where {T <: AbstractFloat}

    !isempty(X) || error("X must be non-empty.")

    vs = Memory{T}(undef, length(partition))
    fill!(vs, zero(T))
    for k in eachindex(partition)

        itr = (X[n] for n in partition[k])
        vs[k] += Statistics.var(itr; corrected = false)
    end
    return vs
end

function compute_partition_for_reduction(
        level_trait::LevelOption,
        dist_callable,
        X::AbstractVector{<:AbstractVector{<:Number}},
    )

    s = SLINKState(get_real_type(X), length(X))
    slink!(s, dist_callable, X)

    level = pick_level(level_trait, s, X)
    return construct_partition(s, level)
end


### frontends

function reduce_pts(
        level_trait::LevelOption,
        center_trait::PartRepOption,
        dist_callable,
        X::AbstractVector{<:AbstractVector{NT}},
    ) where {NT <: Number}

    if length(X) == 1
        return _singleton_reduction(X)
    end
    partition = compute_partition_for_reduction(level_trait, dist_callable, X)

    centers = compute_representatives(center_trait, X, partition)
    variances = _compute_variances(X, partition)

    return centers, variances, partition
end

# the version that also also reduces y.
# usage: reducing training set of scalar-valued (X, y) regression data.
# limit X and y to stride-1, 1-indexing arrays.
function reduce_pts(
        level_trait::LevelOption,
        center_trait::PartRepOption,
        dist_callable,
        X::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}},
        y::AbstractVector{YT},
    ) where {T <: AbstractFloat, YT <: Number}

    if length(X) == 1
        return _singleton_reduction(X)
    end
    partition = compute_partition_for_reduction(level_trait, dist_callable, X)

    Xc = compute_representatives(center_trait, X, partition)
    Xv = _compute_variances(X, partition)

    yc = compute_representatives(center_trait, y, partition)
    yv = _compute_variances(y, partition)

    return Xc, Xv, yc, yv, partition
end

# usage: reducing training set of vector-valued (X, y_set) regression data.
function reduce_pts(
        level_trait::LevelOption,
        center_trait::PartRepOption,
        dist_callable,
        X::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}},
        y_set::AbstractVector{<:AbstractVector{YT}},
    ) where {T <: AbstractFloat, YT <: Number}

    if length(X) == 1
        return _singleton_reduction(X)
    end

    #partition = getreduceptspartition(level_trait, metric, X)
    partition = compute_partition_for_reduction(level_trait, dist_callable, X)

    Xc = compute_representatives(center_trait, X, partition)
    Xv = _compute_variances(X, partition)

    yc_set = collect(compute_representatives(center_trait, y, partition) for y in y_set)
    yv_set = collect(_compute_variances(y, partition) for y in y_set)

    return Xc, Xv, yc_set, yv_set, partition
end

function _singleton_reduction(X::AbstractVector{AbstractVector{T}}) where {T <: Number}
    centers = _create_mem(X)

    variances = Memory{T}(undef, 1)
    fill!(variances, 0)

    partition = Memory{Memory{Int}}(undef, 1)
    partition[begin] = Memory{Int}(undef, 1)
    partition[begin][begin] = 1

    return centers, variances, partition
end

function _create_mem(X::AbstractVector{AbstractVector{T}}) where {T <: Number}

    Z = Memory{Memory{T}}(undef, length(X))
    for i in eachindex(X, Z)
        Z[i] = Memory{T}(X[i])
    end
    return Z
end
