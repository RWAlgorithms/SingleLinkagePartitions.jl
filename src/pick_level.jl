# SPDX-License-Identifier: AGPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

"""
    abstract type LevelOption end

A trait that is a subtype of `LevelOption` is used to pick a level of the partition tree, given a tolerance. Use `subtypes` to see the implemented traits.
"""
abstract type LevelOption end

struct UseCumulativeSLDistance{T} <: LevelOption
    atol::T

    @doc """
        UseCumulativeSLDistance(atol::T) where T <: AbstractFloat

    Creates a variable of type `UseCumulativeSLDistance{T} <: LevelOption`.

    Given a partition tree, the part/cluster that merged as we traverse from level `l` to level `l+1` is the *merging distance* of the two levels.

    `UseCumulativeSLDistance` specifies the threshold on the cumulative merging distance from level `0` to `l`. For use with `reduce_pts`. The internals of `UseCumulativeSLDistance` are considered private API.
    """
    function UseCumulativeSLDistance(atol::T) where {T <: AbstractFloat}
        return new{T}(atol)
    end
end

struct UseSLDistance{T} <: LevelOption
    atol::T

    @doc """
        UseSLDistance(atol::T) where T <: AbstractFloat

    Creates a variable of type `UseSLDistance{T} <: LevelOption`.

    Given a partition tree, the part/cluster that merged as we traverse from level `l` to level `l+1` is the *merging distance* of the two levels.

    `UseSLDistance` specifies the threshold merging distance of level `l` and `l+1`. For use with `reduce_pts`. The internals of `UseSLDistance` are considered private API.
    """
    function UseSLDistance(atol::T) where {T <: AbstractFloat}
        return new{T}(atol)
    end
end

"""
    pick_level(C::UseSLDistance, s::SLINKState, unused_args...)

See `UseSLDistance` for a description of *merging distance*.
Returns the partition for first level `l-1` such that the merging distance between levels ``l`` and ``l-1`` is larger than the `C.atol` struct field. The search starts from ``l`` = 1, ``l-1`` = 0 (the leaf level).

This function allocates. Run `slink!` on `s` first. Return type is `Int`.
"""
function pick_level(C::UseSLDistance, s::SLINKState, unused_args...)

    w = construct_distances(s)

    level = length(w)
    ind = findfirst(xx -> xx > C.atol, w)
    if !isnothing(ind)
        level = ind
    end

    level -= 1
    return clamp(level, 0, get_cardinality(s))
end

"""
    pick_level(C::UseCumulativeSLDistance, s::SLINKState, unused_args...)

See `UseCumulativeSLDistance` for a description of *merging distance*.
Returns the partition for level `l-1` such that the cumulative sum of merging distances from the level `0` (the leaf level) to level `l`is larger than the `C.atol` struct field. The search starts from `l = 1`, `l-1 = 0`.

This function allocates. Run `slink!` on `s` first. Return type is `Int`.
"""
function pick_level(C::UseCumulativeSLDistance, s::SLINKState, unused_args...)

    w = construct_distances(s)

    ind = findfirst(xx -> xx > C.atol, cumsum(w))
    level = length(w)
    if !isnothing(ind)
        level = ind
    end

    level -= 1
    return clamp(level, 0, get_cardinality(s))

    return level
end


######


struct UseMaxDeviation{T <: AbstractFloat} <: LevelOption
    max_dev::T
    max_iters::Int
    search_convergence_atol::T

    @doc """
        UseMaxDeviation(max_dev::T, max_iters::Int, search_convergence_atol::T) where {T <: AbstractFloat}

    Creates a variable of type `UseMaxDeviation{T} <: LevelOption`.

    Terminology:
    Given a point set X on a D-dimensional vector space, we define the *maximum deviation of a partition* `P_X` of `X` as:
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

    If we order the partitions of a partition tree such that the next partition merges two of the parts in the current partition, i.e. in a nesting order, then the sequence of maximum deviations of the ordered partitions is approximately monotonically increasing. This is because a partition with larger sized parts (partitions in the latter parts of the sequence) in general have more variability within their larger sized parts than a partition with smaller sized parts (partitions in the beginning parts of the sequence).

    Description:
    The `UseMaxDeviation` trait/config data type specifies that we should choose a partition from a partition tree such that the chosen partition has a maximum deviation that is approximately the best match to a user-specified target value, which is stored as the `max_dev` struct field.

    The reason the implemented methods associated with `UseMaxDeviation` does not guaranteed to return the best match is because the sequence of maximum deviations that corresponds to the nesting order of partitions is not guaranteed to be monotonic.

    The `max_dev` input should be a non-negative number. It is clamped to zero if the `max_dev` provided is a negative number. Setting this to `0` corresponds to selecting the trivial partition where every point is a singleton part, i.e., level `0` of the partition tree.

    The internals of `UseMaxDeviation` is considered private API.
    """
    function UseMaxDeviation(max_dev::T, max_iters::Int, search_convergence_atol::T) where {T <: AbstractFloat}
        max_dev = clamp(max_dev, zero(T), T(Inf))
        return new{T}(max_dev, max_iters, search_convergence_atol)
    end
end

"""
    UseMaxDeviation(max_dev::T) where {T <: AbstractFloat}

A convenience constructor for `UseMaxDeviation(max_dev, 100, max_dev/100)`.
"""
function UseMaxDeviation(max_dev::T) where {T <: AbstractFloat}
    return UseMaxDeviation(max_dev, 100, max_dev / 100)
end

"""
    pick_level(
        C::UseMaxDeviation{T},
        s::SLINKState{T},
        X::AbstractVector{<:AbstractVector{<:Number}},
    ) where {T}

Returns the level index, of type `Int`, from the partition tree such that its maximum deviation is approximately a best match to the maximum deviation threshold recorded in `C`.

This method uses a bracketing binary search between levels 0 and `get_cardinality(s)-1`. See `UseMaxDeviation` for further details on the terminology.
"""
function pick_level(
        C::UseMaxDeviation{T},
        s::SLINKState{T},
        X::AbstractVector{<:AbstractVector{<:Number}},
    ) where {T}

    max_level = get_cardinality(s) - 1
    if C.max_dev == zero(T)
        return 0
    end

    # binary search.
    f = ll -> compute_max_deviation(X, s, clamp(round(Int, ll), 0, max_level))

    level_float, _ = binary_search(
        f,
        C.max_dev,
        zero(T),
        convert(T, max_level);
        atol = C.search_convergence_atol,
        max_iters = C.max_iters,
    )

    #@show level_float, iters_ran, f(level_float), C.search_convergence_atol

    level = clamp(round(Int, level_float), 0, max_level)
    #@show f(level)
    while f(level) > C.max_dev && level >= 0
        level -= 1
    end

    return level
end

# For this function to return a good approximate to the target_score, `f` must be a strictly increasing monotonic function, and that `f(p) == target_score` for some p such that `p_lb <= p <= p_ub`.
# `f` must be a callable variable that takes an input of type `T` and returns an output of type `T`.
function binary_search(
        f,
        target_score::T,
        p_lb::T,
        p_ub::T;
        atol::T = convert(T, 1.0e-5),
        max_iters = 100,
    ) where {T <: AbstractFloat}

    # initialize.
    lb, ub = p_lb, p_ub
    p = (lb + ub) / 2

    for iter in 1:max_iters
        p = (lb + ub) / 2
        f_p = f(p)

        #@show f_p, target_score, lb, ub
        if isapprox(f_p, target_score; atol = atol) || # found solution.
                isapprox(lb, ub; atol = atol) # we're stuck.

            return p, iter
        end

        if f_p > target_score
            ub = p
        else
            lb = p
        end
    end

    return p, max_iters
end

# computes ρ(X[partition]); see documentation website.
"""
    compute_max_deviation(X, s::SLINKState, level::Integer)

This version extracts the partition assocaited with `level` from `pt`, then calls `compute_max_deviation(X, partition)`.
"""
function compute_max_deviation(X, s::SLINKState, level::Integer)
    return compute_max_deviation(X, construct_partition(s, level))
end

"""
    compute_max_deviation(
        X::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}},
        partition::AbstractVector{<:AbstractVector{Int}},
    ) where T <: AbstractFloat

Given a point set `X` on a *D*-dimensional vector space, we define the *maximum deviation of a partition* `P_X` of `X` as:
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

Computes the maximum of the maximum deviations of a partition of `X`.
This is `ρ_dim(X[partition])` from the documentation website home page.

Returns a number of type `T`.
"""
function compute_max_deviation(
        X::Union{AbstractVector{<:AbstractVector{T}}, AbstractVector{<:AbstractVector{Complex{T}}}},
        partition::AbstractVector{<:AbstractVector{Int}},
    ) where {T <: AbstractFloat}

    !isempty(X) || error("X must be non-empty.")
    D = length(X[begin])

    out = zero(T)
    for k in eachindex(partition)

        μ = Statistics.mean(X[n] for n in partition[k])
        for d in 1:D
            # max norm, aka l-infinity norm of X[n] - μ, for n ∈ partition[k].
            vs_k_d = maximum(abs(X[n][d] - μ[d]) for n in partition[k])
            out = max(vs_k_d, out)
        end
    end
    return out
end
