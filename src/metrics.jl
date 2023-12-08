
abstract type MetricType end

struct DistancesjlMetric{DT} <: MetricType
    dist_metric::DT
end

function geteuclideanmetric()
    #return DistancesjlMetric(Distances.SqEuclidean())
    return DistancesjlMetric(Distances.Euclidean())
end

# assume input is of floating-point type when the metric is from Distances.jl.
# sqrt(dot(x, W*y)) = norm(x-W*y). x, y are real-valued vectors.
struct InnerProductNorm{WT <: Number} <: MetricType
    W::Matrix{WT}
end


struct GeneralMetric{T <: AbstractFloat, MT} <: MetricType
    dist_type::T
    evalmetric::MT # takes two vector inputs, returns a number of type T.
end

function evalmetric(A::GeneralMetric{T,MT}, x::Vector, y::Vector)::T where {T,MT}

    return A.evalmetric(x, y)
end

# XT can be real or complex.
function evalmetric(
    A::DistancesjlMetric,
    x::Union{Vector{Complex{T}}, Vector{T}},
    y::Union{Vector{Complex{T}}, Vector{T}},
    )::T where T <: AbstractFloat
    
    return Distances.evaluate(A.dist_metric, x, y)
end

# XT can be real or complex.
function evalmetric(
    A::InnerProductNorm,
    x::Union{Vector{T}, Vector{Complex{T}}},
    y::Union{Vector{T}, Vector{Complex{T}}},
    )::T where T <: AbstractFloat
    
    return Distances.evaluate(Distances.Euclidean(), x, A.W*y)
end

############# distance computation

# convert to matrix, then use Distances.jl.
function vec2mat(X::Vector{Vector{T}})::Matrix{T} where T
   return hcat(X...)
end

# convinence function only for floating point inputs.
function getpairwisedists(
    metric::MetricType,
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    ) where T <: AbstractFloat

    N = length(X)
    R = Matrix{T}(undef, N, N)
    getpairwisedists!(metric, R, X)

    return R
end

# XT can be real or complex.
function getpairwisedists!(
    metric::GeneralMetric,
    R::Matrix{T}, # mutates, output
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    ) where T <: AbstractFloat

    N = length(X)
    @assert size(R) == (N,N)

    # upper triangle and diagonal.
    for j in axes(R,2)
        for i in Iterators.take(axes(R,1), j)
            R[i,j] = metric.evalmetric(X[i], X[j])
        end
    end

    # lower-triangle.
    for j in axes(R,2)
        for i in Iterators.drop(axes(R,1), j)
            R[i,j] = R[j,i]
        end
    end

    return nothing
end

# For floating point inputs if metric is from Distances.jl.
function getpairwisedists!(
    metric::InnerProductNorm,
    R::Matrix{T}, # mutates, output
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    ) where T <: AbstractFloat

    return getpairwisedists!(metric, R, map(xx->metric.W*xx, X))
end

function getpairwisedists!(
    metric::DistancesjlMetric,
    R::Matrix{T}, # mutates, output
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    ) where T <: AbstractFloat

    return getpairwisedists!(metric, R, vec2mat(X))
end

function getpairwisedists!(
    metric::DistancesjlMetric,
    R::Matrix{T}, # mutates, output
    X_mat::Union{Matrix{T}, Matrix{Complex{T}}},
    ) where T <: AbstractFloat

    Distances.pairwise!(metric.dist_metric, R, X_mat; dims = 2)    
    return nothing
end



############ distance search.

# also used for diagnostics
# returns empty arrays if X is 1-element.
function getdistancesflat(
    metric::MetricType,
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    ) where T <: AbstractFloat

    @assert !isempty(X)
    R = getpairwisedists(metric, X)
    
    return getdistancesflat(R)
end

function getdistancesflat(R::Matrix{T}) where T <: AbstractFloat

    @assert !isempty(R)
    
    N = size(R,1)
    @assert size(R,2) == N

    M = div(N*(N-1), 2)
    A = zeros(T, M)
    inds = Vector{Tuple{Int,Int}}(undef, M)

    k = 0
    for j = 2:N
        for i = 1:j-1 
            k += 1
            #A[k] = norm(X[i]-X[j])
            #A[k] = evalmetric(metric, X[i], X[j])
            A[k] = R[i,j]
            inds[k] = (i,j)
        end
    end

    return A, inds
end