
using LinearAlgebra
import Statistics

using Serialization
using BenchmarkTools

using Revise
import SingleLinkagePartitions
const SL = SingleLinkagePartitions
const Graphs = SL.Graphs
const Distances = SL.Distances