
# assumes f(ps), ps = LinRange(p_lb, p_ub, N) is a monotonic sequence of real numbers for any positive finite integer N.
function binarysearch(
    f::Function,
    target_score::T,
    p_lb::T,
    p_ub::T;
    atol::T = convert(T, 1e-5),
    max_iters = 100,
    ) where T <: AbstractFloat

    # initialize.
    lb, ub = p_lb, p_ub
    p = (lb+ub)/2

    for iter = 1:max_iters
        p = (lb+ub)/2
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

# computes ρ_dim(X[partition]); see the documentation website.
function computedeviationdims(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    partition::Vector{Vector{Int}},
    )::Vector{Vector{T}} where T <: AbstractFloat

    @assert !isempty(X)
    D = length(X[begin])

    vs = Vector{Vector{T}}(undef, length(partition))
    for k in eachindex(partition)

        μ = Statistics.mean( X[n] for n in partition[k] )

        vs[k] = zeros(T, D)
        for d = 1:D
            # max norm, aka l-infinity norm of X[n] - μ, for n ∈ partition[k].
            vs[k][d] = maximum( abs(X[n][d] - μ[d]) for n in partition[k] )
        end
    end

    return vs
end

# computes ρ(X[partition]); see documentation website.
"""
```
computemaxdeviation(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    pt::PartitionTree,
    level::Integer,
) where T <: AbstractFloat
```
This version extracts the partition assocaited with `level` from `pt`, then calls `computemaxdeviation(X, partition)`.
"""
function computemaxdeviation(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    pt::PartitionTree,
    level::Integer,
    ) where T <: AbstractFloat

    return computemaxdeviation(X, getpartition(pt, level))
end

"""
```
computemaxdeviation(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    partition::Vector{Vector{Int}},
) where T <: AbstractFloat
```
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
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    partition::Vector{Vector{Int}},
    ) where T <: AbstractFloat

    ds_X = computedeviationdims(X, partition)
    max_ds = collect( maximum(ds_X[k]) for k in eachindex(ds_X) )
    
    return maximum(max_ds)
end

"""
```
struct UseMaxDeviation{T <: AbstractFloat} <: LevelOption
    max_dev::T # a positive number that specifies a level from a partition tree. For use with `picklevel`.

    # the following are tuning parameters for use with `iteratedsl`.
    max_iters::Int
    search_convergence_atol::T
    verbose::Bool
end
```
Dispatch type/trait, for use with picking a level via `picklevel` or iterated single-linkage via `iteratedsl`.

Terminology:
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

If we order the partitions of a partition tree such that the next partition merges two of the parts in the current partition, i.e. in a nesting order, then the sequence of maximum deviations of the ordered partitions is approximately monotonically increasing. This is because a partition with larger sized parts (partitions in the latter parts of the sequence) in general have more variability within their larger sized parts than a partition with smaller sized parts (partitions in the beginning parts of the sequence).

Description:
The `UseMaxDeviation` trait/config data type specifies that we should choose a partition from a partition tree such that the chosen partition has a maximum deviation that is approximately the best match to a user-specified target value, which is stored as the `max_dev` struct field.

The reason the implemented methods associated with `UseMaxDeviation` does not guaranteed to return the best match is because the sequence of maximum deviations that corresponds to the nesting order of partitions is not guaranteed to be monotonic.
"""
struct UseMaxDeviation{T <: AbstractFloat} <: LevelOption
    max_dev::T # non-negative number. zero means we select the trivial level of 0.
    max_iters::Int
    search_convergence_atol::T
    verbose::Bool
end

"""
```
UseMaxDeviation(max_dev::T)::UseMaxDeviation{T} where T <: AbstractFloat
```
A convenience default constructor.
returns UseMaxDeviation(max_dev, 100, max_dev/100, false)
"""
function UseMaxDeviation(max_dev::T)::UseMaxDeviation{T} where T <: AbstractFloat
    return UseMaxDeviation(max_dev, false)
end

"""
```
UseMaxDeviation(max_dev::T, verbose::Bool)::UseMaxDeviation{T} where T <: AbstractFloat
```
A convenience default constructor.
returns UseMaxDeviation(max_dev, 100, max_dev/100, verbose)
"""
function UseMaxDeviation(max_dev::T, verbose::Bool)::UseMaxDeviation{T} where T <: AbstractFloat
    @assert max_dev >= zero(T) # if zero, picklevel() should return 0.
    # if max_dev is some large number that exceeds maximum(pt.w), then picklevel() should return getNedges(pt)
    return UseMaxDeviation(max_dev, 100, max_dev/100, verbose)
end

"""
```
picklevel(
    C::UseMaxDeviation{T},
    pt::PartitionTree{T},
    X,
)::Int where T <: AbstractFloat
```

Returns an approximate best match partition from `pt` such that its maximum deviation is approximately a best match to `C.max_dev`.
    
Use a bracketing binary search between levels 0 and `getmaxlevel(pt)`. See `UseMaxDeviation` for further details on the terminology.
"""
function picklevel(
    C::UseMaxDeviation{T},
    pt::PartitionTree{T},
    X,
    )::Int where T <: AbstractFloat

    max_level = getNedges(pt)
    if C.max_dev == zero(T)
        return 0
    end

    # if C.max_dev >= pt.w[end] # pt.w should be sorted in ascending order. #maximum(pt.w)
    #     return max_level
    # end

    # binary search.
    f = ll->computemaxdeviation(X, pt, clamp(round(Int,ll), 0, max_level))

    level_float, iters_ran = binarysearch(
        f,
        C.max_dev,
        zero(T),
        convert(T, max_level);
        atol = C.search_convergence_atol,
        max_iters = C.max_iters,
    )
    if C.verbose && iters_ran == C.max_iters
        println("Warning: binarysearch() ran for the specified maximum number of iterations: $iters_ran.")
    end

    #@show level_float, iters_ran, f(level_float), C.search_convergence_atol

    level = clamp(round(Int, level_float), 0, max_level)
    #@show f(level)
    while f(level) > C.max_dev && level >= 0 
        level -= 1
    end

    return level
end

######

function translatepartition!(P::Vector{Vector{Int}}, labels::Vector{Int})
    @assert length( collect(Iterators.flatten(P) )) == length(labels)
    
    for k in eachindex(P)
        for i in eachindex(P[k])
            P[k][i] = labels[P[k][i]]
        end
    end

    return nothing
end

"""
```
iteratedsl(
    level_config::UseMaxDeviation,
    metric::MetricType,
    X0::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}};
    acceptance_factor::T = convert(T, 0.7),
    max_iter = 100,
)::Tuple{Vector{Vector{Int}}, Int} where T <: AbstractFloat
```

Returns a partition `P``, and `iters_ran`.
Iteratively run single-linkage clustering with `UseMaxDeviation ` and `picklevel(` to select a partition at every iteration. Some of the parts in the selected partition are chosen to be moved to a solution partition, and the associated points won't be used for the next iteration of single-linkage clustering.
See the documentation website for details on the iteration algorithm, including terminating conditions.

"""
function iteratedsl(
    level_config::UseMaxDeviation,
    metric::MetricType,
    X0::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}};
    acceptance_factor::T = convert(T, 0.7),
    max_iter = 100,
    )::Tuple{Vector{Vector{Int}}, Int} where T <: AbstractFloat

    @assert zero(T) < acceptance_factor < one(T)

    #max_dev = level_config.max_dev
    #keep_dev = max_dev*acceptance_factor

    X = X0
    X_labels = collect(1:length(X0))

    P = Vector{Vector{Int}}(undef, 0) # output partition.

    for iters = 1:max_iter

        if isempty(X)
            return P, iters
        end

        if length(X) == 1
            push!(P, [X_labels[begin];] )
            return P, iters
        end

        # compute partition tree.
        pt = computesl(metric, X)
        level = picklevel(level_config, pt, X)
        partition = getpartition(pt, level)
        

        ds_X = computedeviationdims(X, partition)
        max_ds = collect( maximum(ds_X[k]) for k in eachindex(ds_X) )
        max_all_ds = maximum(max_ds)
        #max_all_ds = computemaxdeviation(X, partition)
        keep_dev = acceptance_factor*max_all_ds

        # these inds are for partition.
        update_inds = findall(xx->xx>=keep_dev, max_ds) # LHS can't be empty due to previous if statement.
        next_inds = setdiff(1:length(partition), update_inds)

        #@show length(update_inds), keep_dev, max_all_ds # debug.

        if isempty(next_inds)
            # no points left.

            # update output partition, P.
            translatepartition!(partition, X_labels)
            append!(P, partition[update_inds])

            return P, iters
        end

        # this is for the pt set, X.
        next_X_inds = collect(
            Iterators.flatten(partition[next_inds])
        )

        # update output partition, P.
        translatepartition!(partition, X_labels)
        append!(P, partition[update_inds])
        #@show P

        # prepare next iteration's input, X.
        X = X[next_X_inds]
        X_labels = X_labels[next_X_inds]
    end

    # last round.
    pt = computesl(metric, X)
    level = picklevel(level_config, pt, X)

    partition = getpartition(pt, level)
    translatepartition!(partition, X_labels)
    append!(P, partition)

    return P, max_iter
end