
# this is a minimum-spanning tree. Edges sorted in ascending order.
struct PartitionTree{T <: AbstractFloat}
    
    N_nodes::Int

    # edges.
    src_nodes::Vector{Int}
    dest_nodes::Vector{Int}
    w::Vector{T}
end

function getdistances(pt::PartitionTree{T})::Vector{T} where T <: AbstractFloat
    return pt.w
end

function generateallpartitions(pt::PartitionTree)
    return collect(
        getpartition(pt, level)
        for level = 0:getNedges(pt)
    )
end

function getNedges(pt::PartitionTree)::Int
    return length(pt.w)
end

function checkdatastructure(pt::PartitionTree)::Bool

    flag1 = pt.N_nodes > 0
    flag2 = length(pt.w) == length(pt.dest_nodes) == length(pt.src_nodes)
    flag3 = length(pt.w) > 0
    flag4 = issorted(pt.w)

    return flag1 && flag2 && flag3 && flag4
end

function PartitionTree(R::Matrix{T}) where T <: AbstractFloat

    @assert size(R,1) == size(R,2)
    
    mst, w_mst = computemst(R)

    inds = sortperm(w_mst)
    src_nodes = collect(Graphs.src(mst[k]) for k in inds)
    dest_nodes = collect(Graphs.dst(mst[k]) for k in inds)

    return PartitionTree(size(R,1), src_nodes, dest_nodes, w_mst[inds])
end

# this is the front-end to single-linkage partition computation.
function computesl(
    metric::Union{DistancesjlMetric, InnerProductNorm},
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    )::PartitionTree{T} where T <: AbstractFloat

    R = getpairwisedists(metric, X)
    # R_flat, R_flat_inds = SL.getdistancesflat(R)
    # #sorted_R = sort(R_flat)

    return PartitionTree(R)
end

function computemst(metric::MetricType, X::Vector{Vector{T}}) where T

    R = getpairwisedists(metric, X)
    return computemst(R)
end

function computemst(R::Matrix{T}) where T <: AbstractFloat

    N = size(R,1)
    @assert size(R,2) == N

    g = Graphs.complete_graph(N)

    mst = Graphs.prim_mst(g, R)

    # # construct minimum spanning tree.
    w_mst = collect( 
        R[Graphs.src(mst[k]), Graphs.dst(mst[k])]
        for k in eachindex(mst)
    )

    return mst, w_mst #, g, tree
end

function getgraph(pt::PartitionTree)
    
    @assert checkdatastructure(pt)

    tree_graph = Graphs.SimpleGraph(pt.N_nodes)
    for k in eachindex(pt.src_nodes)
        Graphs.add_edge!(tree_graph, pt.src_nodes[k], pt.dest_nodes[k])
    end

    return tree_graph
end

# for diagnostics.
# does not care about the connectivity of the minimal-spanning tree, mst.
function findedge(mst, src::Integer, dest::Integer)
    return findfirst(
        xx->(
            (Graphs.src(xx)==src && Graphs.dst(xx)==dest) ||
            (Graphs.dst(xx)==src && Graphs.src(xx)==dest)
        ),
        mst,
    )
end

# level 0 is the leaf level. level == length(mst) is the root level, where all nodes are in the singlton part.
# level indicates the number of merges between two subgraphs, from the elaf level.
#  the leaf level is 0.
function getpartition(pt::PartitionTree, level::Integer)::Vector{Vector{Int}}

    max_level = getNedges(pt) # the root level.

    @assert 0 <= level <= max_level
    @assert checkdatastructure(pt)

    # manually prepare the leaf and root levels to save computation.
    if level == 0
        # Leaf level: this partition has all singleton parts.
        return collect( [i;] for i = 1:pt.N_nodes)
    end
    
    if level == max_level
        # Root level: this partition has only one part, which has all pts.
        partition = Vector{Vector{Int}}(undef, 1)
        partition[begin] = collect(1:pt.N_nodes)

        return partition
    end

    # The other levels.
    g = Graphs.SimpleGraph(pt.N_nodes)
    for k in Iterators.take(eachindex(pt.src_nodes), level)
        Graphs.add_edge!(g, pt.src_nodes[k], pt.dest_nodes[k])
    end

    partition = Graphs.connected_components(g)
    return partition
end

