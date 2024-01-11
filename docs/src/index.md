# SingleLinkagePartitions.jl

This package implements the single-linkage clustering algorithm via a minimum spanning tree approach. It is suitable for handling data that can fit in the local machine's RAM.

Single-linkage clustering is a hierarchical clustering method, which means for a dataset of ``N`` points in ``\mathbb{R}^{D}``, it generates a set of ``N`` nested partitions. The implementation here uses a minimum-spanning tree, which should offer some efficiency and flexibility than the naive single-linkage clustering algorithm.

The documentation is on-going, but demo walk-throughs for the chaining example as well as the quick-start examples in the README.MD file should address most basic usage cases. A walkthrough for merging points might be authored in the future.


# Terminology and algorithms

I devised some strategies to select a partition given a computed single-linkage partition tree. This section describes the terminologies and algorithms associated with these strategies.

## Pick a partition based on maximum deviation

Let the dimension of the input array be ``D``. Given a positive integer ``L``, the ``L``-th order Riesz transform is a collection of iterated Riesz transforms, each of which is specified by a multi-index ``a \in \mathbb{A}\left(D,L\right)``. The index set is defined as:

The **dimensional maximum deviation from the mean of a point set** ``A\subset\mathbb{R}^{D}`` is defined as: 
```math
\begin{align*}
\rho_{\text{dim}}\left(A\right) & :=\left\{ \max_{a\in A}\left\{ \left|a_{d}-\mu_{A,d}\right|\right\} \right\} _{d\in\left[D\right]}\cong\begin{bmatrix}\max_{a\in A}\left\{ \left|a_{1}-\mu_{A,1}\right|\right\} \\
\max_{a\in A}\left\{ \left|a_{2}-\mu_{A,2}\right|\right\} \\
\vdots\\
\max_{a\in A}\left\{ \left|a_{D}-\mu_{A,D}\right|\right\} 
\end{bmatrix},\\
\bar{\rho}_{\text{dim}}\left(A\right) & =\max\rho_{\text{dim}}\left(A\right)\\
 & =\max_{d\in\left[D\right],\,a\in A}\left\{ \left|a_{d}-\mu_{A,d}\right|\right\} .
\end{align*}
```
The subscripts mean the coordinates of the point ``a\in A``, and ``\mu_{A}`` is the average of the points in ``A``.

Let ``X\subset\mathbb{R}^{D}`` be a finite point set, and ``P_{X}`` be a partition of ``X``. The **maximum deviation from the mean of a partition** ``P_{X}`` of a point set ``X`` is defined as:
```math
\begin{align*}
\rho\left(P_{X}\right) & :=\max_{U\in P_{X}}\left\{ \bar{\rho}_{\text{dim}}\left(U\right)\right\} \\
 & =\max_{U\in P_{X},d\in\left[D\right],a\in U}\left|a_{d}-\mu_{U,d}\right|.
\end{align*}
```
Each partition tree level ``l\in\left\{ 0,1,\cdots,L\right\}`` corresponds to a partition, which we denote by ``P_{X,l}``. The squence ``\left\{ \rho\left(P_{X,l}\right)\right\} _{l}`` is approximately monotonically increasing, so a bracketed binary search algorithm is used to find a partition that approximately best matches a specified target ``\rho_{\text{target}}``.

After the binary search, if the selected partition ``X_{l}`` has tree index ``l``, and exceeds ``\rho_{\text{target}}``, then we the sequence ``k\in\left\{ l-1,l-2,\cdots,0\right\} `` until we have ``\rho\left(P_{X_{l}}\right)<\rho_{\text{target}}``. Therefore, we gaurantee that the select partition has a maximum deviation that is less than ``\rho_{\text{target}}``.



### Iteration scheme
When we want to reduce the chaining behavior, we could iteratively restart single-linkage clustering for smaller and smaller subsets of the original point set ``X``. My algorithm is as follows:

- Initialize ``R_{0}`` to the empty set, and ``X_{0}:=X``.
- For each iteration ``j``, we want to compute the ``j``-th iteration
update set ``M_{j}``, the results set of subsets ``R_{j}``, and the
process set ``X_{j}``. 

For the ``j``-th iteration:
1. Pick a partition ``P_{X_{j}}`` from the tree via the maximum deviation strategy discussed earlier.
2. Define the threshold ``\tau_{j}=\beta\cdot\rho\left(P_{X_{j}}\right)``, where the user-specific **acceptance factor** ``\beta\in\left(0,1\right)`` is used to make sure ``\tau_{j}`` approaches a maximum of ``\rho\left(P_{X_{j}}\right)``.
3. Perform the update:
```math
\begin{align*}
M_{j} & :=\left\{ U\in P_{X_{j}}\,\mid\,\bar{\rho}_{\text{dim}}\left(U\right)\geq\tau_{j}\right\} ,\\
R_{j} & :=R_{j-1}\cup M_{j},\\
X_{j} & :=X_{j-1}\setminus M_{j}.
\end{align*}
```
Check for the following stopping conditions. Increment ``j`` and repeat from step 1 if conditions are not met. Upon exit of this iterative procedure, ``R_{j}`` for the termination iteration ``j`` should be a partition of ``X``.

- If ``X_{j}`` is empty, then exit with no post-processing.

- If ``X_{j}`` has only one point, then add that point as a singleton subset into ``R_{j}``, then exit.

- If ``M_{j}`` is empty, then assign ``R_{j}=R_{j-1}\cup P_{X_{j}}`` and exit.

# Citation
You can use the *Cite This Repository* button below the *About* section on the GitHub webpage of this repository.