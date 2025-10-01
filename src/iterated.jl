# SPDX-License-Identifier: AGPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

function translate_partition!(P::AbstractVector{<:AbstractVector{Int}}, labels::AbstractVector{Int})
    @assert length(collect(Iterators.flatten(P))) == length(labels)

    for k in eachindex(P)
        for i in eachindex(P[k])
            P[k][i] = labels[P[k][i]]
        end
    end

    return nothing
end

"""
    iterated_sl(
        level_config::UseMaxDeviation,
        dist_callable,
        X0::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}};
        acceptance_factor::T = convert(T, 0.7),
        max_iter::Int = 100,
    ) where {T <: AbstractFloat}


Returns a partition `P``, and `iters_ran`. Return type is `Tuple{Vector{Vector{Int}}, Int}`.
Iteratively run single-linkage clustering with `UseMaxDeviation ` and `picklevel(` to select a partition at every iteration. Some of the parts in the selected partition are chosen to be moved to a solution partition, and the associated points won't be used for the next iteration of single-linkage clustering.
See the documentation website for details on the iteration algorithm, including terminating conditions.

"""
function iterated_sl(
        level_config::UseMaxDeviation,
        dist_callable,
        X0::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}};
        acceptance_factor::T = convert(T, 0.7),
        max_iter::Int = 100,
    ) where {T <: AbstractFloat}

    @assert zero(T) < acceptance_factor < one(T)

    X = X0
    X_labels = collect(1:length(X0))

    P = Vector{Vector{Int}}(undef, 0) # output partition.

    for iters in 1:max_iter

        if isempty(X)
            return P, iters
        end

        if length(X) == 1
            push!(P, [X_labels[begin];])
            return P, iters
        end

        # compute partition tree.
        # pt = computesl(metric, X)
        # level = picklevel(level_config, pt, X)
        # partition = getpartition(pt, level)
        partition = compute_partition_for_reduction(level_config, dist_callable, X)

        ds_X = compute_deviation_dims(X, partition)
        max_ds = collect(maximum(ds_X[k]) for k in eachindex(ds_X))
        max_all_ds = maximum(max_ds)
        #max_all_ds = computemaxdeviation(X, partition)
        keep_dev = acceptance_factor * max_all_ds

        # these inds are for partition.
        update_inds = findall(xx -> xx >= keep_dev, max_ds) # LHS can't be empty due to previous if statement.
        next_inds = setdiff(1:length(partition), update_inds)

        if isempty(next_inds)
            # no points left.

            # update output partition, P.
            translate_partition!(partition, X_labels)
            append!(P, partition[update_inds])

            return P, iters
        end

        # this is for the pt set, X.
        next_X_inds = collect(
            Iterators.flatten(partition[next_inds])
        )

        # update output partition, P.
        translate_partition!(partition, X_labels)
        append!(P, partition[update_inds])
        #@show P

        # prepare next iteration's input, X.
        X = X[next_X_inds]
        X_labels = X_labels[next_X_inds]
    end

    # last round.
    # pt = computesl(metric, X)
    # level = picklevel(level_config, pt, X)
    # partition = getpartition(pt, level)
    partition = compute_partition_for_reduction(level_config, dist_callable, X)

    translate_partition!(partition, X_labels)
    append!(P, partition)

    return P, max_iter
end

# computes ρ_dim(X[partition]); see the documentation website.
# return type: Vector{Vector{T}}
function compute_deviation_dims(
        X::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}},
        partition::AbstractVector{<:AbstractVector{Int}},
    ) where {T <: AbstractFloat}

    @assert !isempty(X)
    D = length(X[begin])

    vs = Vector{Vector{T}}(undef, length(partition))
    for k in eachindex(partition)

        μ = Statistics.mean(X[n] for n in partition[k])

        vs[k] = zeros(T, D)
        for d in 1:D
            # max norm, aka l-infinity norm of X[n] - μ, for n ∈ partition[k].
            vs[k][d] = maximum(abs(X[n][d] - μ[d]) for n in partition[k])
        end
    end

    return vs
end


# computes ρ(X[partition]); see documentation website.
"""
    computemaxdeviation(
        X::Union{AbstractVector{<: AbstractVector{T}}, AbstractVector{ <: AbstractVector{Complex{T}}}},
        s::SLINKState,
        level::Integer,
    ) where T <: AbstractFloat

This version extracts the partition assocaited with `level` from `pt`, then calls `computemaxdeviation(X, partition)`.
"""
function computemaxdeviation(
        X::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}},
        s::SLINKState,
        level::Integer,
    ) where {T <: AbstractFloat}

    return compute_deviation_dims(X, getpartition(pt, level))
end


"""
    computemaxdeviation(
        X::Union{AbstractVector{<: AbstractVector{T}}, AbstractVector{ <: AbstractVector{Complex{T}}}},
        partition::Vector{Vector{Int}},
    ) where T <: AbstractFloat

Given a point set X on a D-dimensional vector space, we define the *maximum deviation of a partition* P_X of X as:
```julia
ρ = maximum(
    maximum(
        maximum(
            abs(a[d]-mean(A)[d])
            for a in A
        )
        for A in P_X
    )
    for d in 1:D
)
```

`computemaxdeviation` computes ρ a bit more efficiently, with
'''julia
P_X = collect( X[partition[k]] for k in eachindex(partition) )
'''
"""
function computemaxdeviation(
        X::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}},
        partition::Vector{Vector{Int}},
    ) where {T <: AbstractFloat}

    ds_X = computedeviationdims(X, partition)
    max_ds = collect(maximum(ds_X[k]) for k in eachindex(ds_X))

    return maximum(max_ds)
end
