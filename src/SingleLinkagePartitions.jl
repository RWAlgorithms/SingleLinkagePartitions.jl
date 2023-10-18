module SingleLinkagePartitions

#include("utils.jl")  # TODO add convinence method so that built-in distance metrics like the ones in utils.jl can be used. Also allow users to wrap function in a custom struct type for custom emtricfuncs.

include("core.jl")
#include("merge.jl") # TODO clean up the interface for mergepoints() and impleemnt efficint algorithms for merge points.

include("front_end.jl")

export runsinglelinkage,
    mergepoints,
    mergepointsfull,
    getdistances,
    instantiatepartition,
    fuseparts

end # module SingleLinkagePartitions
