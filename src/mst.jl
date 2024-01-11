
# this is a minimum-spanning tree. Edges sorted in ascending order.
struct PartitionTree{T <: AbstractFloat}
    
    N_nodes::Int

    # edges.
    src_nodes::Vector{Int}
    dest_nodes::Vector{Int}
    w::Vector{T} # merging distances.
end

"""
```
getdistances(pt::PartitionTree{T})::Vector{T} where T <: AbstractFloat
```
Returns the set of minimum distanes between consecutive nested partitions.
If `pt` was constructed from `N` points, then there are `N-1` entries in the returned 1-D array of distances.
"""
function getdistances(pt::PartitionTree{T})::Vector{T} where T <: AbstractFloat
    return pt.w
end

"""
```
getlargestdistance(pt::PartitionTree{T})::Vector{T} where T <: AbstractFloat
```
Returns the largest merge/single-linkage distance.
"""
function getlargestdistance(pt::PartitionTree{T})::Vector{T} where T <: AbstractFloat
    return pt.w[end] # sorted in ascending order.
end

"""
```
generateallpartitions(pt::PartitionTree)::Vector{Vector{Vector{Int}}}
```
returns the set of all partitions in the partition tree, `partition_set`
If `N` points was used to generate the partition tree `pt`, then:
- `partition_set` has `N+1` entries. Each entry is a partition.
- Given the input point set that generated `pt`, call it `X::Vector{Vector{T}}`, a partition index 1 <= m <= N+1, part index k, part member index i, if we assign `z = partition_set[m][k][i]`, then `X[z]` is point that belongs to the m-th partition, k-th part/cluster in the partition, and it is the i-th point in that part/cluster.

"""
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

"""
```
computesl(
    metric::Union{DistancesjlMetric, InnerProductNorm},
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    )::PartitionTree{T} where T <: AbstractFloat
```
Computes and returns the single-linkage partition tree of `X`, with respect to distance metric, `metric`.
The partitions in the tree are specified by a *level* index from 0 to `maxlevel(pt)`, where `pt` is the returned `PartitionTree` variable.

See `getpartition` to select a partition by its level from the return partition tree.
"""
function computesl(
    metric::Union{DistancesjlMetric, InnerProductNorm},
    X::Union{Vector{Vector{T}}, Vector{Vector{Complex{T}}}},
    )::PartitionTree{T} where T <: AbstractFloat

    R = getpairwisedists(metric, X)

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

"""
```
getpartition(pt::PartitionTree, level::Integer)::Vector{Vector{Int}}
```

`level` indicates the number of merges that occured to produce a partition. Think of it as an index from 0 to `max_level = getmaxlevel(pt)`.

We use the terminology *part* to mean an element (i.e. a subset) of a partition. Some folks call it a *cluster*.

`level` == 0 corresponds to the leaf level where each point is assigned its own part (i.e. cluster), i.e., we have only singleton parts.

The maximum level, `level` == `getmaxlevel(pt)`, corresponds to the root level where all points are assigned to be in the same part.
"""
function getpartition(pt::PartitionTree, level::Integer)::Vector{Vector{Int}}

    max_level = getNedges(pt) # the root level.

    @assert 0 <= level <= max_level
    X_is_non_singleton = checkdatastructure(pt)
    if !X_is_non_singleton
        return collect( ones(Int, 1) for _ = 1:1 )
    end

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
