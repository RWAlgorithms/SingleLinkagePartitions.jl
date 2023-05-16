module SingleLinkagePartitions

include("core.jl")
include("front_end.jl")

export runsinglelinkage,
    mergepoints,
    mergepointsfull,
    getdistances,
    instantiatepartition,
    fuseparts

end # module SingleLinkagePartitions
