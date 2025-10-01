# SPDX-License-Identifier: AGPL-3.0-only
# Copyright © 2025 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

"""
    get_max_level(s::SLINKState)

Returns the maximum level of the partition tree that `s` represents. Run `slink!` on `s` first.
"""
function get_max_level(s::SLINKState)
    return get_cardinality(s) - 1
end

"""
    avg_duplicates(X::AbstractVector, y::AbstractVector, zero_tol::AbstractFloat)

For any locations in X that are nearly identical (with respect to `zero_tol`), reduce the corresponding `y` entries by `mean`. Use `replace_duplicates` if the points to be merged should be determined by a score-based approach.

`zero_tol` is the single-linkage distance threshold. It is similar to a tolerance for merging the points.

Returns
- `Xc`: a reduced version of `X`.
- `yc`: a reduced version of `y`.
"""
function avg_duplicates(X::AbstractVector, y::AbstractVector, zero_tol::AbstractFloat)
    Xc, _, yc, _ = reduce_pts(
        UseSLDistance(zero_tol),
        UseMean(),
        EuclideanDistance(),
        X,
        y,
    )
    return Xc, yc
end

"""
    replace_duplicates(
        copy_trait::CopyTrait,
        score_trait::ScoreTrait,
        X::AbstractVector,
        y::AbstractVector{<:AbstractFloat},
        zero_tol::AbstractFloat,
    )

For any locations in X that are nearly identical (with respect to `zero_tol`), reduce the corresponding `y` entries by the behavior specified by `score_trait`. Use `avg_duplicates` if the mean of the points to be merged should be used instead of a score-based approach.

`zero_tol` is the single-linkage distance threshold. It is similar to a tolerance for merging the points.

Returns
- `Xc`: a reduced version of `X`.
- `yc`: a reduced version of `y`.
"""
function replace_duplicates(
        copy_trait::CopyTrait,
        score_trait::ScoreTrait,
        X::AbstractVector,
        y::AbstractVector{<:AbstractFloat},
        zero_tol::AbstractFloat,
    )

    # when Xc can be any pt from which X is constructed from.
    # score_trait specifies the type of replacement for the duplicates.
    center_trait = UseScore(copy_trait, score_trait, y)
    level_trait = UseSLDistance(zero_tol)

    Xc, _, yc, _ = reduce_pts(
        level_trait,
        center_trait,
        EuclideanDistance(),
        X,
        y,
    )
    return Xc, yc
end
