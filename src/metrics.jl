
abstract type MetricType end

# a specific metric vs.itself squared doesn't actually matter for singlelinkage. Use this instead of Euclidean.
struct EuclideanSquared <: MetricType end

struct Euclidean <: MetricType end

struct WeightedEuclidean{T} <: MetricType
    G::Matrix{T}
end

struct LpNorm{T} <: MetricType
    p::T
end

struct GeneralMetric{MT} <: MetricType
    evalmetric::MT # takes two vector inputs, returns a real number.
end

function evalmetric(::Euclidean, x::Vector{T}, y::Vector{T})::T where T

    out = evalmetric(x, y, metric)

    return forcepositive(sqrt(out))
end

function evalmetric(::EuclideanSquared, x::Vector{T}, y::Vector{T})::T where T

    out = zero(T)
    for i in eachindex(x)
        out += (x[i] - y[i])^2
    end

    return out
end

function evalmetric(::LpNorm, x::Vector{T}, y::Vector{T})::T where T

    return norm(x-y, p) # TODO devectorize.
end

function evalmetric(A::WeightedEuclidean, x::Vector{T}, y::Vector{T})::T where T

    r = x-y
    out = dot(r, A.G*r) # TODO devectorize this.

    return forcepositive(sqrt(out))
end

function evalmetric(A::GeneralMetric, x::Vector{T}, y::Vector{T})::T where T

    return A.evalmetric(x, y)
end