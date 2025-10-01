# SPDX-License-Identifier: AGPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

module SingleLinkagePartitions

using LinearAlgebra
using Statistics

"""
    vecs_to_mat(X::AbstractVector{<:AbstractVector{T}}) where {T}

Returns a matrix constructed from `X`, with `X[d,n] == X[n][d]`. Return type is `Matrix{T}`.
"""
function vecs_to_mat(X::AbstractVector{<:AbstractVector{T}}) where {T}
    !isempty(X) || error("X must be non-empty.")
    return reshape(collect(Iterators.flatten(X)), length(X[begin]), length(X))
end

function get_real_type(::AbstractVector{<:AbstractVector{Complex{T}}}) where {T <: AbstractFloat}
    return T
end

function get_real_type(::AbstractVector{<:AbstractVector{T}}) where {T <: AbstractFloat}
    return T
end

function get_real_type(::AbstractMatrix{T}) where {T <: AbstractFloat}
    return T
end

function get_num_entries(A::AbstractMatrix)
    return size(A, 2)
end

function get_num_entries(x::AbstractVector)
    return length(x)
end

include("slink.jl")
include("pick_level.jl")
include("reduce.jl")
include("frontend.jl")

include("iterated.jl")

public slink!, EuclideanDistance, SLINKState,
    construct_distances, construct_partition, construct_partition_tree,
    get_cardinality, vecs_to_mat,

    pick_level, get_max_level,

    LevelOption,
    UseCumulativeSLDistance, UseSLDistance, UseMaxDeviation,

    PartRepOption,
    UseProximitytoMean, UseMean, # UseScore.

    UseScore,
    UseMaximum, UseMinimum, ScoreTrait,

    reduce_pts, avg_duplicates, replace_duplicates,

    compute_max_deviation, iterated_sl

end # module SingleLinkagePartitions
