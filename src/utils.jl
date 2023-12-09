
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

# maximum diviations from the mean.
function getmaxdeviations(
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
            vs[k][d] = maximum( abs(X[n][d] - μ[d]) for n in partition[k] )
        end
    end

    return vs
end

function computemaxdeviation(
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    pt::PartitionTree,
    level::Integer,
    ) where T <: AbstractFloat

    partition = getpartition(pt, level)
    ds_X = getmaxdeviations(X, partition)
    max_ds = collect( maximum(ds_X[k]) for k in eachindex(ds_X) )
    
    return maximum(max_ds)
end


# The pts in a part cannot have a magnitude deviation from its mean that is larger than max_dev.
struct UseMaxDeviation{T <: AbstractFloat} <: LevelOption
    max_dev::T # non-negative number. zero means we select the trivial level of 0.
    max_iters::Int
    search_convergence_atol::T
    verbose::Bool
end

function UseMaxDeviation(max_dev::T)::UseMaxDeviation{T} where T <: AbstractFloat
    return UseMaxDeviation(max_dev, false)
end

function UseMaxDeviation(max_dev::T, verbose::Bool)::UseMaxDeviation{T} where T <: AbstractFloat
    @assert max_dev >= zero(T) # if zero, getlevel() should return 0.
    # if max_dev is some large number that exceeds maximum(pt.w), then getlevel() should return getNedges(pt)
    return UseMaxDeviation(max_dev, 100, max_dev/100, verbose)
end


function getlevel(
    C::UseMaxDeviation{T},
    pt::PartitionTree{T},
    X,
    )::Int where T <: AbstractFloat
    

    max_level = getNedges(pt)
    if C.max_dev == zero(T)
        return 0
    end

    if C.max_dev >= pt.w[end] # pt.w should be sorted in ascending order. #maximum(pt.w)
        return max_level
    end

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
        println("Warning: binarysearch() ran for max iters.")
    end

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


function iteratedsl(
    level_config::UseMaxDeviation,
    metric::MetricType,
    X0::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}};
    acceptance_factor::T = convert(T, 0.7),
    max_iter = 100,
    ) where T <: AbstractFloat

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
        level = getlevel(level_config, pt, X)
        partition = getpartition(pt, level)
        

        ds_X = getmaxdeviations(X, partition)
        max_ds = collect( maximum(ds_X[k]) for k in eachindex(ds_X) )

        max_all_ds = maximum(max_ds)
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
        @show P

        # prepare next iteration's input, X.
        X = X[next_X_inds]
        X_labels = X_labels[next_X_inds]
    end

    # last round.
    pt = computesl(metric, X)
    level = getlevel(level_config, pt, X)

    partition = getpartition(pt, level)
    translatepartition!(partition, X_labels)
    append!(P, partition)

    return P, max_iter
end