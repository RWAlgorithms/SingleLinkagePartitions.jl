module SingleLinkagePartitions


using LinearAlgebra
using Statistics

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

#include("utils.jl")  # TODO add convinence method so that built-in distance metrics like the ones in utils.jl can be used. Also allow users to wrap function in a custom struct type for custom emtricfuncs.

include("core.jl")
include("consolidate.jl")
# TODO clean up the interface for mergepoints() and impleemnt efficint algorithms for merge points.

include("front_end.jl")

export runsinglelinkage,
    mergepoints,
    mergepointsfull,
    getdistances,
    instantiatepartition,
    fuseparts

end # module SingleLinkagePartitions
