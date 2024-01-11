module SingleLinkagePartitions


using LinearAlgebra
using Statistics
import Graphs
import Distances

function forcepositive(x::Float32)::Float32
    return clamp(x, 0.0f0, Inf32)
end

function forcepositive(x::Float64)::Float64
    return clamp(x, 0.0, Inf)
end

function forcepositive(x::Type{T})::T where T <: AbstractFloat
    return clamp(x, zero(T), convert(T, Inf))
end

include("metrics.jl")
include("mst.jl")

include("reduce.jl")
include("iterated.jl")

include("frontend.jl")

export getmaxlevel, #

    UseImportance,
    UseProximitytoMean,
    UseMean,

    #DistancesjlMetric,
    #InnerProductNorm,
    #GeneralMetric,
    geteuclideanmetric, #
    UseMaxDeviation, #
    UseSLDistance, #
    UseCumulativeSLDistance,

    picklevel, #
    computemaxdeviation, #
    iteratedsl, #

    reducepts,
    computesl, #
    generateallpartitions, #
    
    #getpairwisedists,
    #getpairwisedists!,
    getpartition, #
    getdistances, #
    getlargestdistance #
    
end # module SingleLinkagePartitions
