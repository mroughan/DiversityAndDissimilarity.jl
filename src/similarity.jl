"""Abstract supertype for pairwise similarity, dissimilarity, and distance indices."""
abstract type PairwiseIndex <: DiversityIndex end

"""
Jaccard incidence similarity.

```math
J(A,B) = \\frac{|A \\cap B|}{|A \\cup B|}
```
"""
struct Jaccard <: PairwiseIndex end

"""
Sorensen-Dice incidence similarity.

```math
S(A,B) = \\frac{2|A \\cap B|}{|A| + |B|}
```
"""
struct SorensenDice <: PairwiseIndex end

"""
Overlap (Szymkiewicz-Simpson) incidence similarity.

```math
O(A,B) = \\frac{|A \\cap B|}{\\min(|A|, |B|)}
```
"""
struct Overlap <: PairwiseIndex end

"""
Bray-Curtis abundance dissimilarity.

```math
BC(x,y) = \\frac{\\sum_i |x_i - y_i|}{\\sum_i (x_i + y_i)}
```
"""
struct BrayCurtis <: PairwiseIndex end

"""
Ruzicka, or quantitative Jaccard, abundance similarity.

```math
R(x,y) = \\frac{\\sum_i \\min(x_i,y_i)}{\\sum_i \\max(x_i,y_i)}
```
"""
struct Ruzicka <: PairwiseIndex end

"""
Total variation distance between probability vectors.

```math
TV(p,q) = \\frac{1}{2}\\sum_i |p_i - q_i|
```
"""
struct TotalVariation <: PairwiseIndex end

"""
Manhattan, or L1, distance between probability vectors.

```math
d_1(p,q) = \\sum_i |p_i - q_i|
```
"""
struct Manhattan <: PairwiseIndex end

"""
Euclidean, or L2, distance between probability vectors.

```math
d_2(p,q) = \\sqrt{\\sum_i (p_i - q_i)^2}
```
"""
struct Euclidean <: PairwiseIndex end

"""
Averaged Canberra distance between abundance vectors.

```math
C(x,y) = \\frac{1}{m}\\sum_{i:x_i+y_i>0}\\frac{|x_i-y_i|}{x_i+y_i}
```
"""
struct Canberra <: PairwiseIndex end

"""
Hellinger distance between probability vectors.

```math
H(p,q) = \\frac{1}{\\sqrt{2}}\\sqrt{\\sum_i (\\sqrt{p_i} - \\sqrt{q_i})^2}
```
"""
struct Hellinger <: PairwiseIndex end

"""
Chord distance between square-root transformed probability vectors.

```math
d_c(p,q) = \\sqrt{\\sum_i (\\sqrt{p_i} - \\sqrt{q_i})^2}
```
"""
struct Chord <: PairwiseIndex end

"""
Bhattacharyya coefficient and distance between probability vectors.

```math
BC(p,q) = \\sum_i \\sqrt{p_i q_i}, \\qquad
d_B(p,q) = -\\log BC(p,q)
```
"""
struct Bhattacharyya <: PairwiseIndex end

"""
Kullback-Leibler divergence between probability vectors.

`dissimilarity(KullbackLeibler(), left, right)` returns ``D_{KL}(p \\Vert q)``:

```math
D_{KL}(p \\Vert q) = \\sum_i p_i \\log_b \\frac{p_i}{q_i}
```

**This divergence is asymmetric**: `dissimilarity(KullbackLeibler(), a, b)` and
`dissimilarity(KullbackLeibler(), b, a)` are generally different, and a full
community distance matrix will not be symmetric. See [`is_symmetric`](@ref).

Use `estimator` for low-sample corrections. Supported options mirror Shannon
entropy estimation: [`MillerMadow`](@ref), [`AddGamma`](@ref) for pseudocounts
(`AddGamma(1)` is Laplace and `AddGamma(0.5)` is Jeffreys), [`HausserStrimmer`](@ref)
for shrinkage, and [`ChaoShen`](@ref) for a Good-Turing unseen-mass correction.
Without smoothing or unseen-mass correction this divergence returns `Inf` when
`p_i > 0` and `q_i == 0` for any coordinate.

The [`MillerMadow`](@ref) correction subtracts the standard entropy bias
correction ``(S-1)/(2n \\log b)`` from the plugin divergence estimate, where
``S`` is the number of observed positive-probability categories in `left` and
``n`` is its total count. This is a first-order correction for the bias in
``H(p)``; the cross-entropy term ``H(p,q)`` is left uncorrected.
"""
struct KullbackLeibler{E<:ShannonEstimator,S} <: PairwiseIndex
    base::Float64
    estimator::E
    support::S
end

function KullbackLeibler(; base::Real=2, estimator::ShannonEstimator=Plugin(), support=nothing)
    base > 0 && base != 1 || throw(ArgumentError("base must be positive and not equal to 1"))
    return KullbackLeibler(float(base), estimator, support)
end

"""
Absolute Shannon entropy difference between probability vectors.

```math
|H_b(p) - H_b(q)|
```

!!! warning "Does not measure distributional divergence"
    This index compares the scalar entropy magnitudes of two assemblages; it
    does not measure how different the distributions themselves are. Two
    assemblages with completely disjoint species but identical species-abundance
    profiles will score zero. Use [`JensenShannon`](@ref) or
    [`KullbackLeibler`](@ref) when you want a proper distributional divergence.
"""
struct ShannonDifference{E<:ShannonEstimator,S} <: PairwiseIndex
    base::Float64
    estimator::E
    support::S
end

function ShannonDifference(; base::Real=2, estimator::ShannonEstimator=Plugin(), support=nothing)
    base > 0 && base != 1 || throw(ArgumentError("base must be positive and not equal to 1"))
    return ShannonDifference(float(base), estimator, support)
end

"""
Jensen difference of Shannon entropy between probability vectors.

```math
J_H(p,q) = H_b\\left(\\frac{p+q}{2}\\right) - \\frac{H_b(p)+H_b(q)}{2}
```

For Shannon entropy this equals the Jensen-Shannon divergence. It returns the
raw divergence value. Use [`JensenShannon`](@ref) when you want the metric
square-root form (i.e. `dissimilarity(JensenShannon(), ...)` with `distance=true`).
Use `estimator` for the same low-sample corrections available to [`KullbackLeibler`](@ref).
"""
struct JensenDifference{E<:ShannonEstimator,S} <: PairwiseIndex
    base::Float64
    estimator::E
    support::S
end

function JensenDifference(; base::Real=2, estimator::ShannonEstimator=Plugin(), support=nothing)
    base > 0 && base != 1 || throw(ArgumentError("base must be positive and not equal to 1"))
    return JensenDifference(float(base), estimator, support)
end

"""
Jensen-Shannon divergence or distance between probability vectors.

`JensenShannon(; base=2, distance=true)` returns the square root of the
Jensen-Shannon divergence, which is a proper metric. Set `distance=false` for
the divergence itself (identical to [`JensenDifference`](@ref)). Use `estimator`
for low-sample corrections such as `MillerMadow()`, `AddGamma(1)` for Laplace
smoothing, `AddGamma(0.5)` for Jeffreys smoothing, `HausserStrimmer()` shrinkage,
and `ChaoShen()` / Good-Turing unseen-mass correction.

See also [`JensenDifference`](@ref), which returns the raw divergence value
without taking a square root.
"""
struct JensenShannon{E<:ShannonEstimator,S} <: PairwiseIndex
    base::Float64
    distance::Bool
    estimator::E
    support::S
end

function JensenShannon(; base::Real=2, distance::Bool=true,
        estimator::ShannonEstimator=Plugin(), support=nothing)
    base > 0 && base != 1 || throw(ArgumentError("base must be positive and not equal to 1"))
    return JensenShannon(float(base), distance, estimator, support)
end

"""
Morisita-Horn abundance similarity.

```math
MH(x,y) =
\\frac{2\\sum_i x_i y_i}{(\\lambda_x + \\lambda_y) N_x N_y},
\\qquad
\\lambda_x = \\frac{\\sum_i x_i^2}{N_x^2}
```
"""
struct MorisitaHorn <: PairwiseIndex end

"""
    similarity(index, left, right; frequencies=true)

Compare two assemblages using an incidence or abundance similarity index.
"""
function similarity(::Jaccard, left, right; frequencies::Bool=true)
    left_species = _species_set(left; frequencies)
    right_species = _species_set(right; frequencies)
    denominator = length(union(left_species, right_species))
    denominator == 0 && throw(ArgumentError("at least one species is required"))
    return length(intersect(left_species, right_species)) / denominator
end

function similarity(::SorensenDice, left, right; frequencies::Bool=true)
    left_species = _species_set(left; frequencies)
    right_species = _species_set(right; frequencies)
    denominator = length(left_species) + length(right_species)
    denominator == 0 && throw(ArgumentError("at least one species is required"))
    return 2length(intersect(left_species, right_species)) / denominator
end

function similarity(::Overlap, left, right; frequencies::Bool=true)
    left_species = _species_set(left; frequencies)
    right_species = _species_set(right; frequencies)
    denominator = min(length(left_species), length(right_species))
    denominator == 0 && throw(ArgumentError("both assemblages must contain at least one species"))
    return length(intersect(left_species, right_species)) / denominator
end

function dissimilarity(::BrayCurtis, left, right; frequencies::Bool=true)
    left_abundance, right_abundance = _checked_aligned_abundances(left, right; frequencies)
    denominator = sum(left_abundance) + sum(right_abundance)
    denominator > 0 || throw(ArgumentError("total abundance must be positive"))
    return sum(abs, left_abundance .- right_abundance) / denominator
end

function dissimilarity(::BrayCurtis, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing)
    _check_matrix_frequencies(frequencies)
    _validate_community_matrix(data)
    nsites = size(data, 1)
    ntaxa = size(data, 2)
    abundance = permutedims(data)
    row_totals = Vector{Float64}(undef, nsites)
    for row in 1:nsites
        total = 0.0
        @inbounds for column in 1:ntaxa
            total += abundance[column, row]
        end
        row_totals[row] = total
    end

    result = Matrix{Float64}(undef, nsites, nsites)
    for i in 1:nsites
        result[i, i] = 0.0
        for j in (i + 1):nsites
            numerator = 0.0
            @inbounds for column in 1:ntaxa
                numerator += abs(abundance[column, i] - abundance[column, j])
            end
            value = numerator / (row_totals[i] + row_totals[j])
            @inbounds begin
                result[i, j] = value
                result[j, i] = value
            end
        end
    end
    return result
end

function similarity(::Ruzicka, left, right; frequencies::Bool=true)
    left_abundance, right_abundance = _checked_aligned_abundances(left, right; frequencies)
    denominator = sum(max.(left_abundance, right_abundance))
    denominator > 0 || throw(ArgumentError("total abundance must be positive"))
    return sum(min.(left_abundance, right_abundance)) / denominator
end

function dissimilarity(::TotalVariation, left, right; frequencies::Bool=true)
    left_probability, right_probability = _aligned_probabilities(left, right; frequencies)
    return 0.5 * sum(abs, left_probability .- right_probability)
end

function dissimilarity(::Manhattan, left, right; frequencies::Bool=true)
    left_probability, right_probability = _aligned_probabilities(left, right; frequencies)
    return sum(abs, left_probability .- right_probability)
end

function dissimilarity(::Euclidean, left, right; frequencies::Bool=true)
    left_probability, right_probability = _aligned_probabilities(left, right; frequencies)
    return sqrt(sum(abs2, left_probability .- right_probability))
end

function dissimilarity(::Canberra, left, right; frequencies::Bool=true)
    left_abundance, right_abundance = _checked_aligned_abundances(left, right; frequencies)
    total = 0.0
    terms = 0
    for (left_value, right_value) in zip(left_abundance, right_abundance)
        denominator = left_value + right_value
        if denominator > 0
            total += abs(left_value - right_value) / denominator
            terms += 1
        end
    end
    terms > 0 || throw(ArgumentError("total abundance must be positive"))
    return total / terms
end

function dissimilarity(::Hellinger, left, right; frequencies::Bool=true)
    left_probability, right_probability = _aligned_probabilities(left, right; frequencies)
    return sqrt(sum(abs2, sqrt.(left_probability) .- sqrt.(right_probability)) / 2)
end

function dissimilarity(::Chord, left, right; frequencies::Bool=true)
    left_probability, right_probability = _aligned_probabilities(left, right; frequencies)
    return sqrt(sum(abs2, sqrt.(left_probability) .- sqrt.(right_probability)))
end

function similarity(::Bhattacharyya, left, right; frequencies::Bool=true)
    left_probability, right_probability = _aligned_probabilities(left, right; frequencies)
    return sum(sqrt.(left_probability .* right_probability))
end

function dissimilarity(index::Bhattacharyya, left, right; frequencies::Bool=true)
    coefficient = similarity(index, left, right; frequencies)
    return coefficient == 0 ? Inf : -log(coefficient)
end

function dissimilarity(index::KullbackLeibler, left, right; frequencies::Bool=true)
    left_counts, right_counts, left_total, right_total =
        _aligned_divergence_counts(left, right; frequencies, support=index.support)
    left_probability, right_probability = _divergence_probability_pair(
        index.estimator, left_counts, right_counts, left_total, right_total)
    divergence = _kl_divergence(
        left_probability,
        right_probability,
        index.base,
    )
    if index.estimator isa MillerMadow
        correction = (count(>(0), left_counts) - 1) / (2left_total * log(index.base))
        divergence = max(0.0, divergence - correction)
    end
    return divergence
end

function dissimilarity(index::KullbackLeibler, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing)
    _check_matrix_frequencies(frequencies)
    _validate_community_matrix(data)
    nsites = size(data, 1)
    result = Matrix{Float64}(undef, nsites, nsites)
    for i in 1:nsites
        result[i, i] = 0.0
        for j in 1:nsites
            if i != j
                result[i, j] = dissimilarity(index, view(data, i, :), view(data, j, :); frequencies)
            end
        end
    end
    return result
end

function dissimilarity(index::KullbackLeibler, data; frequencies::Bool=true, species=nothing)
    matrix = community_matrix(data; species)
    return dissimilarity(index, matrix; frequencies)
end

function dissimilarity(index::ShannonDifference, left, right; frequencies::Bool=true)
    left_counts, right_counts, left_total, right_total =
        _aligned_divergence_counts(left, right; frequencies, support=index.support)
    return abs(
        _divergence_entropy(index.estimator, left_counts, left_total, index.base) -
        _divergence_entropy(index.estimator, right_counts, right_total, index.base),
    )
end

function similarity(index::ShannonDifference, left, right; frequencies::Bool=true)
    left_counts, right_counts, _, _ =
        _aligned_divergence_counts(left, right; frequencies, support=index.support)
    support = _effective_divergence_support(index.estimator, left_counts, right_counts)
    support <= 1 && return 1.0
    max_difference = log(support) / log(index.base)
    return 1 - dissimilarity(index, left, right; frequencies) / max_difference
end

function _jensen_difference(index::JensenDifference, left, right; frequencies::Bool=true)
    left_counts, right_counts, left_total, right_total =
        _aligned_divergence_counts(left, right; frequencies, support=index.support)
    return _jensen_difference(index.estimator, left_counts, right_counts,
        left_total, right_total, index.base)
end

function dissimilarity(index::JensenDifference, left, right; frequencies::Bool=true)
    return _jensen_difference(index, left, right; frequencies)
end

function similarity(index::JensenDifference, left, right; frequencies::Bool=true)
    max_difference = log(2) / log(index.base)
    return 1 - dissimilarity(index, left, right; frequencies) / max_difference
end

function _jensen_shannon_divergence(index::JensenShannon, left, right; frequencies::Bool=true)
    left_counts, right_counts, left_total, right_total =
        _aligned_divergence_counts(left, right; frequencies, support=index.support)
    return _jensen_difference(index.estimator, left_counts, right_counts,
        left_total, right_total, index.base)
end

function dissimilarity(index::JensenShannon, left, right; frequencies::Bool=true)
    divergence = _jensen_shannon_divergence(index, left, right; frequencies)
    return index.distance ? sqrt(divergence) : divergence
end

function similarity(index::JensenShannon, left, right; frequencies::Bool=true)
    max_distance = index.distance ? sqrt(log(2) / log(index.base)) : log(2) / log(index.base)
    return 1 - dissimilarity(index, left, right; frequencies) / max_distance
end

function similarity(::MorisitaHorn, left, right; frequencies::Bool=true)
    left_abundance, right_abundance = _checked_aligned_abundances(left, right; frequencies)
    left_total = sum(left_abundance)
    right_total = sum(right_abundance)
    left_total > 0 && right_total > 0 || throw(ArgumentError("both assemblages must have positive total abundance"))
    left_lambda = sum(abs2, left_abundance) / left_total^2
    right_lambda = sum(abs2, right_abundance) / right_total^2
    denominator = (left_lambda + right_lambda) * left_total * right_total
    denominator > 0 || throw(ArgumentError("Morisita-Horn denominator must be positive"))
    return 2sum(left_abundance .* right_abundance) / denominator
end

"""
    dissimilarity(index, left, right; frequencies=true)

Compare two assemblages using a dissimilarity or distance form of an index.
"""
dissimilarity(index::Jaccard, left, right; frequencies::Bool=true) =
    1 - similarity(index, left, right; frequencies)
dissimilarity(index::SorensenDice, left, right; frequencies::Bool=true) =
    1 - similarity(index, left, right; frequencies)
dissimilarity(index::Overlap, left, right; frequencies::Bool=true) =
    1 - similarity(index, left, right; frequencies)
similarity(index::BrayCurtis, left, right; frequencies::Bool=true) =
    1 - dissimilarity(index, left, right; frequencies)
dissimilarity(index::Ruzicka, left, right; frequencies::Bool=true) =
    1 - similarity(index, left, right; frequencies)
similarity(index::TotalVariation, left, right; frequencies::Bool=true) =
    1 - dissimilarity(index, left, right; frequencies)
similarity(index::Canberra, left, right; frequencies::Bool=true) =
    1 - dissimilarity(index, left, right; frequencies)
similarity(index::Hellinger, left, right; frequencies::Bool=true) =
    1 - dissimilarity(index, left, right; frequencies)
dissimilarity(index::MorisitaHorn, left, right; frequencies::Bool=true) =
    1 - similarity(index, left, right; frequencies)

function similarity(::Jaccard, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing)
    _check_matrix_frequencies(frequencies)
    _validate_community_matrix(data)
    presence = _incidence_bitsets(data)
    nsites = size(presence, 2)
    nwords = size(presence, 1)
    result = Matrix{Float64}(undef, nsites, nsites)
    for i in 1:nsites
        result[i, i] = 1.0
        for j in (i + 1):nsites
            intersection = 0
            union_count = 0
            @inbounds for word in 1:nwords
                left = presence[word, i]
                right = presence[word, j]
                intersection += count_ones(left & right)
                union_count += count_ones(left | right)
            end
            value = intersection / union_count
            @inbounds begin
                result[i, j] = value
                result[j, i] = value
            end
        end
    end
    return result
end

function dissimilarity(::Jaccard, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing)
    _check_matrix_frequencies(frequencies)
    _validate_community_matrix(data)
    presence = _incidence_bitsets(data)
    nsites = size(presence, 2)
    nwords = size(presence, 1)
    result = Matrix{Float64}(undef, nsites, nsites)
    for i in 1:nsites
        result[i, i] = 0.0
        for j in (i + 1):nsites
            intersection = 0
            union_count = 0
            @inbounds for word in 1:nwords
                left = presence[word, i]
                right = presence[word, j]
                intersection += count_ones(left & right)
                union_count += count_ones(left | right)
            end
            value = 1 - intersection / union_count
            @inbounds begin
                result[i, j] = value
                result[j, i] = value
            end
        end
    end
    return result
end

function _incidence_bitsets(data::AbstractMatrix{<:Real})
    nsites = size(data, 1)
    ntaxa = size(data, 2)
    nwords = cld(ntaxa, 64)
    presence = zeros(UInt64, nwords, nsites)
    @inbounds for column in 1:ntaxa
        word = ((column - 1) >>> 6) + 1
        bit = UInt64(1) << ((column - 1) & 63)
        for row in 1:nsites
            if data[row, column] > 0
                presence[word, row] |= bit
            end
        end
    end
    return presence
end

function dissimilarity(::Hellinger, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing)
    _check_matrix_frequencies(frequencies)
    _validate_community_matrix(data)
    nsites = size(data, 1)
    ntaxa = size(data, 2)
    abundance = permutedims(data)
    roots = Matrix{Float64}(undef, ntaxa, nsites)
    for row in 1:nsites
        total = 0.0
        @inbounds for column in 1:ntaxa
            total += abundance[column, row]
        end
        @inbounds for column in 1:ntaxa
            roots[column, row] = sqrt(abundance[column, row] / total)
        end
    end

    result = Matrix{Float64}(undef, nsites, nsites)
    for i in 1:nsites
        result[i, i] = 0.0
        for j in (i + 1):nsites
            total = 0.0
            @inbounds for column in 1:ntaxa
                total += abs2(roots[column, i] - roots[column, j])
            end
            value = sqrt(total / 2)
            @inbounds begin
                result[i, j] = value
                result[j, i] = value
            end
        end
    end
    return result
end

"""
    similarity(index, data; frequencies=true, species=nothing)

Return a pairwise similarity matrix across rows of a community matrix or
Tables.jl-compatible table.
"""
function similarity(index::DiversityIndex, data; frequencies::Bool=true, species=nothing)
    matrix = community_matrix(data; species)
    return _pairwise_matrix((left, right) -> similarity(index, left, right; frequencies), matrix)
end

"""
    distance(index, left, right; frequencies=true)

Alias for [`dissimilarity`](@ref).
"""
distance(index, left, right; frequencies::Bool=true) =
    dissimilarity(index, left, right; frequencies)

"""
    dissimilarity(index, data; frequencies=true, species=nothing)

Return a pairwise dissimilarity matrix across rows of a community matrix or
Tables.jl-compatible table.
"""
function dissimilarity(index::DiversityIndex, data; frequencies::Bool=true, species=nothing)
    matrix = community_matrix(data; species)
    return _pairwise_matrix((left, right) -> dissimilarity(index, left, right; frequencies), matrix)
end

"""
    distance(index, data; frequencies=true, species=nothing)

Return a pairwise distance/dissimilarity matrix across rows of a community
matrix or Tables.jl-compatible table.
"""
distance(index::DiversityIndex, data; frequencies::Bool=true, species=nothing) =
    dissimilarity(index, data; frequencies, species)

"""
    labeled_distance(index, data; labels=nothing, label=nothing, frequencies=true, species=nothing)

Return pairwise distances with sample labels as a named tuple
`(labels=..., matrix=...)`.

Pass `labels` explicitly for matrices. For Tables.jl-compatible inputs, pass
`label` as the column name containing sample/site identifiers.
"""
function labeled_distance(index::DiversityIndex, data; labels=nothing, label=nothing,
        frequencies::Bool=true, species=nothing)
    return (
        labels=_pairwise_labels(data; labels, label),
        matrix=distance(index, data; frequencies, species),
    )
end

"""
    labeled_dissimilarity(index, data; labels=nothing, label=nothing, frequencies=true, species=nothing)

Return pairwise dissimilarities with sample labels.
"""
function labeled_dissimilarity(index::DiversityIndex, data; labels=nothing, label=nothing,
        frequencies::Bool=true, species=nothing)
    return (
        labels=_pairwise_labels(data; labels, label),
        matrix=dissimilarity(index, data; frequencies, species),
    )
end

"""
    labeled_similarity(index, data; labels=nothing, label=nothing, frequencies=true, species=nothing)

Return pairwise similarities with sample labels.
"""
function labeled_similarity(index::DiversityIndex, data; labels=nothing, label=nothing,
        frequencies::Bool=true, species=nothing)
    return (
        labels=_pairwise_labels(data; labels, label),
        matrix=similarity(index, data; frequencies, species),
    )
end

function _pairwise_labels(data; labels=nothing, label=nothing)
    labels !== nothing && label !== nothing &&
        throw(ArgumentError("pass either labels or label, not both"))
    if labels !== nothing
        result = collect(labels)
    elseif label !== nothing
        _is_table(data) || throw(ArgumentError("label column selection requires a Tables.jl-compatible input"))
        columns = Tables.columns(data)
        names = collect(Tables.columnnames(columns))
        label in names || throw(ArgumentError("label column $(repr(label)) was not found"))
        result = collect(Tables.getcolumn(columns, label))
    else
        result = collect(1:_sample_count(data))
    end
    length(result) == _sample_count(data) ||
        throw(ArgumentError("number of labels must match the number of samples"))
    return result
end

_sample_count(data::AbstractMatrix) = size(data, 1)
function _sample_count(data)
    if _is_table(data)
        columns = Tables.columns(data)
        names = collect(Tables.columnnames(columns))
        isempty(names) && throw(ArgumentError("table must have at least one column"))
        return length(Tables.getcolumn(columns, first(names)))
    end
    return size(community_matrix(data), 1)
end

function _pairwise_matrix(metric, matrix::AbstractMatrix{<:Real})
    _check_matrix_frequencies(true)
    _validate_community_matrix(matrix)
    n = size(matrix, 1)
    result = Matrix{Float64}(undef, n, n)
    for i in 1:n
        result[i, i] = metric(view(matrix, i, :), view(matrix, i, :))
        for j in (i + 1):n
            value = metric(view(matrix, i, :), view(matrix, j, :))
            result[i, j] = value
            result[j, i] = value
        end
    end
    return result
end

function _checked_aligned_abundances(left, right; frequencies::Bool=true)
    left_abundance, right_abundance = _aligned_abundances(left, right; frequencies)
    length(left_abundance) == length(right_abundance) || throw(DimensionMismatch("abundance vectors must have the same length"))
    for value in Iterators.flatten((left_abundance, right_abundance))
        value < 0 && throw(ArgumentError("abundances must be non-negative"))
        isfinite(value) || throw(ArgumentError("abundances must be finite"))
    end
    sum(left_abundance) > 0 && sum(right_abundance) > 0 ||
        throw(ArgumentError("both assemblages must have positive total abundance"))
    return left_abundance, right_abundance
end

function _aligned_probabilities(left, right; frequencies::Bool=true)
    left_abundance, right_abundance = _checked_aligned_abundances(left, right; frequencies)
    return left_abundance ./ sum(left_abundance), right_abundance ./ sum(right_abundance)
end

function _aligned_divergence_counts(left, right; frequencies::Bool=true, support=nothing)
    if support === nothing
        left_counts, right_counts = _checked_aligned_abundances(left, right; frequencies)
    elseif support isa Integer
        left_counts, right_counts = _checked_aligned_abundances(left, right; frequencies)
        support >= length(left_counts) || throw(ArgumentError("support size must be at least the aligned observed support size"))
        append!(left_counts, zeros(Float64, support - length(left_counts)))
        append!(right_counts, zeros(Float64, support - length(right_counts)))
    else
        known = collect(support)
        left_observed = _divergence_observed_counts(left; frequencies)
        right_observed = _divergence_observed_counts(right; frequencies)
        observed_set = union(Set(keys(left_observed)), Set(keys(right_observed)))
        known_set = Set(known)
        observed_set ⊆ known_set || throw(ArgumentError("support must include every observed category"))
        left_counts = [float(get(left_observed, category, 0)) for category in known]
        right_counts = [float(get(right_observed, category, 0)) for category in known]
        any(v -> v < 0 || !isfinite(v), Iterators.flatten((left_counts, right_counts))) &&
            throw(ArgumentError("abundances must be non-negative and finite"))
    end
    left_total = sum(left_counts)
    right_total = sum(right_counts)
    left_total > 0 && right_total > 0 || throw(ArgumentError("both assemblages must have positive total abundance"))
    return left_counts, right_counts, left_total, right_total
end

function _divergence_observed_counts(data; frequencies::Bool=true)
    if data isa AbstractVector{<:Real} && frequencies
        for value in data
            value < 0 && throw(ArgumentError("abundances must be non-negative"))
            isfinite(value) || throw(ArgumentError("abundances must be finite"))
        end
        return Dict(i => float(value) for (i, value) in pairs(data))
    end
    return counts(data)
end

function _divergence_probabilities(::Union{Plugin,MillerMadow}, counts, n)
    return counts ./ n
end

function _divergence_probabilities(estimator::AddGamma, counts, n)
    support_size = length(counts)
    denominator = n + estimator.gamma * support_size
    return [(count + estimator.gamma) / denominator for count in counts]
end

function _divergence_probabilities(::HausserStrimmer, counts, n)
    support_size = length(counts)
    target = inv(support_size)
    probabilities = counts ./ n
    denominator = sum(p -> (target - p)^2, probabilities)
    shrinkage = if n <= 1 || denominator <= eps(Float64)
        1.0
    else
        numerator = sum(p -> p * (1 - p), probabilities) / (n - 1)
        clamp(numerator / denominator, 0.0, 1.0)
    end
    return [shrinkage * target + (1 - shrinkage) * p for p in probabilities]
end

function _divergence_probabilities(::ChaoShen, counts, n)
    singletons = count(==(1), counts)
    coverage = 1 - singletons / n
    coverage >= 0 || throw(ArgumentError("Good-Turing coverage estimate must be non-negative"))
    observed = [coverage * count / n for count in counts]
    push!(observed, max(0.0, 1 - coverage))
    return observed
end

function _divergence_entropy(estimator::ShannonEstimator, counts, n, base)
    probabilities = _divergence_probabilities(estimator, counts, n)
    return _shannon_entropy(probabilities, base)
end

function _divergence_probability_pair(estimator::ShannonEstimator, left_counts, right_counts,
        left_total, right_total)
    return (
        _divergence_probabilities(estimator, left_counts, left_total),
        _divergence_probabilities(estimator, right_counts, right_total),
    )
end

function _divergence_probability_pair(::ChaoShen, left_counts, right_counts,
        left_total, right_total)
    return (
        _good_turing_pair_probabilities(left_counts, right_counts, left_total),
        _good_turing_pair_probabilities(right_counts, left_counts, right_total),
    )
end

function _good_turing_pair_probabilities(counts, other_counts, n)
    singletons = count(==(1), counts)
    coverage = 1 - singletons / n
    coverage >= 0 || throw(ArgumentError("Good-Turing coverage estimate must be non-negative"))
    unseen_mass = max(0.0, 1 - coverage)
    zero_in_this = [i for i in eachindex(counts) if counts[i] == 0 && other_counts[i] > 0]
    unseen_share = unseen_mass / (length(zero_in_this) + 1)
    probabilities = [count > 0 ? coverage * count / n : 0.0 for count in counts]
    for index in zero_in_this
        probabilities[index] = unseen_share
    end
    push!(probabilities, unseen_share)
    return probabilities
end

function _jensen_difference(estimator::ShannonEstimator, left_counts, right_counts,
        left_total, right_total, base)
    left_probability, right_probability =
        _divergence_probability_pair(estimator, left_counts, right_counts, left_total, right_total)
    mixture = (left_probability .+ right_probability) ./ 2
    value = _shannon_entropy(mixture, base) -
        (_shannon_entropy(left_probability, base) + _shannon_entropy(right_probability, base)) / 2
    if estimator isa MillerMadow
        mixture_total = left_total + right_total
        mixture_support = count(>(0), mixture)
        left_support = count(>(0), left_counts)
        right_support = count(>(0), right_counts)
        correction = (mixture_support - 1) / (2mixture_total * log(base)) -
            0.5 * (
                (left_support - 1) / (2left_total * log(base)) +
                (right_support - 1) / (2right_total * log(base))
            )
        value += correction
    end
    return max(0.0, value)
end

function _effective_divergence_support(estimator::ShannonEstimator, left_counts, right_counts)
    if estimator isa ChaoShen
        return count(value -> value != 0, left_counts .+ right_counts) + 1
    end
    return length(left_counts)
end

function _kl_divergence(left_probability, right_probability, base)
    total = 0.0
    log_base = log(base)
    for (left_value, right_value) in zip(left_probability, right_probability)
        if left_value > 0
            right_value > 0 || return Inf
            total += left_value * log(left_value / right_value) / log_base
        end
    end
    return total
end

function _shannon_entropy(probability, base)
    total = 0.0
    log_base = log(base)
    for value in probability
        if value > 0
            total -= value * log(value) / log_base
        end
    end
    return total
end

"""
    jaccard_index(left, right; frequencies=true)

Return Jaccard similarity between two assemblages.
"""
jaccard_index(left, right; frequencies::Bool=true) = similarity(Jaccard(), left, right; frequencies)
jaccard_index(data; frequencies::Bool=true, species=nothing) =
    similarity(Jaccard(), data; frequencies, species)

"""
    jaccard_similarity(left, right; frequencies=true)

Alias for [`jaccard_index`](@ref).
"""
jaccard_similarity(left, right; frequencies::Bool=true) = jaccard_index(left, right; frequencies)
jaccard_similarity(data; frequencies::Bool=true, species=nothing) =
    jaccard_index(data; frequencies, species)

"""
    jaccard_distance(left, right; frequencies=true)

Return Jaccard dissimilarity between two assemblages.
"""
jaccard_distance(left, right; frequencies::Bool=true) = dissimilarity(Jaccard(), left, right; frequencies)
jaccard_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Jaccard(), data; frequencies, species)

"""
    sorensen_dice_index(left, right; frequencies=true)

Return Sorensen-Dice similarity between two assemblages.
"""
sorensen_dice_index(left, right; frequencies::Bool=true) = similarity(SorensenDice(), left, right; frequencies)
sorensen_dice_index(data; frequencies::Bool=true, species=nothing) =
    similarity(SorensenDice(), data; frequencies, species)

"""
    sorensen_index(left, right; frequencies=true)

Alias for [`sorensen_dice_index`](@ref).
"""
sorensen_index(left, right; frequencies::Bool=true) = sorensen_dice_index(left, right; frequencies)
sorensen_index(data; frequencies::Bool=true, species=nothing) =
    sorensen_dice_index(data; frequencies, species)

"""
    sorensen_dice_dissimilarity(left, right; frequencies=true)

Return Sorensen-Dice dissimilarity between two assemblages.
"""
sorensen_dice_dissimilarity(left, right; frequencies::Bool=true) =
    dissimilarity(SorensenDice(), left, right; frequencies)
sorensen_dice_dissimilarity(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(SorensenDice(), data; frequencies, species)

"""
    sorensen_dice_distance(left, right; frequencies=true)

Alias for [`sorensen_dice_dissimilarity`](@ref).
"""
sorensen_dice_distance(left, right; frequencies::Bool=true) =
    sorensen_dice_dissimilarity(left, right; frequencies)
sorensen_dice_distance(data; frequencies::Bool=true, species=nothing) =
    sorensen_dice_dissimilarity(data; frequencies, species)

"""
    sorensen_distance(left, right; frequencies=true)

Alias for [`sorensen_dice_distance`](@ref).
"""
sorensen_distance(left, right; frequencies::Bool=true) =
    sorensen_dice_distance(left, right; frequencies)
sorensen_distance(data; frequencies::Bool=true, species=nothing) =
    sorensen_dice_distance(data; frequencies, species)

"""
    bray_curtis_dissimilarity(left, right; frequencies=true)

Return Bray-Curtis dissimilarity between two assemblages.
"""
bray_curtis_dissimilarity(left, right; frequencies::Bool=true) =
    dissimilarity(BrayCurtis(), left, right; frequencies)
bray_curtis_dissimilarity(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(BrayCurtis(), data; frequencies, species)

"""
    bray_curtis_distance(left, right; frequencies=true)

Alias for [`bray_curtis_dissimilarity`](@ref).
"""
bray_curtis_distance(left, right; frequencies::Bool=true) =
    bray_curtis_dissimilarity(left, right; frequencies)
bray_curtis_distance(data; frequencies::Bool=true, species=nothing) =
    bray_curtis_dissimilarity(data; frequencies, species)

"""
    overlap_similarity(left, right; frequencies=true)

Return overlap (Szymkiewicz-Simpson) similarity between two assemblages.
"""
overlap_similarity(left, right; frequencies::Bool=true) = similarity(Overlap(), left, right; frequencies)
overlap_similarity(data; frequencies::Bool=true, species=nothing) =
    similarity(Overlap(), data; frequencies, species)

"""
    overlap_distance(left, right; frequencies=true)

Return one minus overlap similarity.
"""
overlap_distance(left, right; frequencies::Bool=true) = dissimilarity(Overlap(), left, right; frequencies)
overlap_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Overlap(), data; frequencies, species)

"""
    ruzicka_similarity(left, right; frequencies=true)

Return Ruzicka, or quantitative Jaccard, abundance similarity.
"""
ruzicka_similarity(left, right; frequencies::Bool=true) = similarity(Ruzicka(), left, right; frequencies)
ruzicka_similarity(data; frequencies::Bool=true, species=nothing) =
    similarity(Ruzicka(), data; frequencies, species)

"""
    quantitative_jaccard_similarity(left, right; frequencies=true)

Alias for [`ruzicka_similarity`](@ref).
"""
quantitative_jaccard_similarity(left, right; frequencies::Bool=true) =
    ruzicka_similarity(left, right; frequencies)
quantitative_jaccard_similarity(data; frequencies::Bool=true, species=nothing) =
    ruzicka_similarity(data; frequencies, species)

"""
    ruzicka_distance(left, right; frequencies=true)

Return one minus Ruzicka similarity.
"""
ruzicka_distance(left, right; frequencies::Bool=true) = dissimilarity(Ruzicka(), left, right; frequencies)
ruzicka_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Ruzicka(), data; frequencies, species)

"""
    quantitative_jaccard_distance(left, right; frequencies=true)

Alias for [`ruzicka_distance`](@ref).
"""
quantitative_jaccard_distance(left, right; frequencies::Bool=true) =
    ruzicka_distance(left, right; frequencies)
quantitative_jaccard_distance(data; frequencies::Bool=true, species=nothing) =
    ruzicka_distance(data; frequencies, species)

"""
    total_variation_distance(left, right; frequencies=true)

Return total variation distance between normalized abundance/probability vectors.
"""
total_variation_distance(left, right; frequencies::Bool=true) =
    dissimilarity(TotalVariation(), left, right; frequencies)
total_variation_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(TotalVariation(), data; frequencies, species)

"""
    manhattan_distance(left, right; frequencies=true)

Return Manhattan/L1 distance between normalized abundance/probability vectors.
"""
manhattan_distance(left, right; frequencies::Bool=true) =
    dissimilarity(Manhattan(), left, right; frequencies)
manhattan_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Manhattan(), data; frequencies, species)

"""
    euclidean_distance(left, right; frequencies=true)

Return Euclidean/L2 distance between normalized abundance/probability vectors.
"""
euclidean_distance(left, right; frequencies::Bool=true) =
    dissimilarity(Euclidean(), left, right; frequencies)
euclidean_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Euclidean(), data; frequencies, species)

"""
    canberra_distance(left, right; frequencies=true)

Return averaged Canberra distance between abundance vectors.
"""
canberra_distance(left, right; frequencies::Bool=true) = dissimilarity(Canberra(), left, right; frequencies)
canberra_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Canberra(), data; frequencies, species)

"""
    hellinger_distance(left, right; frequencies=true)

Return Hellinger distance between normalized abundance/probability vectors.
"""
hellinger_distance(left, right; frequencies::Bool=true) =
    dissimilarity(Hellinger(), left, right; frequencies)
hellinger_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Hellinger(), data; frequencies, species)

"""
    chord_distance(left, right; frequencies=true)

Return chord distance between square-root transformed probabilities.
"""
chord_distance(left, right; frequencies::Bool=true) = dissimilarity(Chord(), left, right; frequencies)
chord_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Chord(), data; frequencies, species)

"""
    bhattacharyya_coefficient(left, right; frequencies=true)

Return the Bhattacharyya coefficient between probability vectors.
"""
bhattacharyya_coefficient(left, right; frequencies::Bool=true) =
    similarity(Bhattacharyya(), left, right; frequencies)
bhattacharyya_coefficient(data; frequencies::Bool=true, species=nothing) =
    similarity(Bhattacharyya(), data; frequencies, species)

"""
    bhattacharyya_distance(left, right; frequencies=true)

Return Bhattacharyya distance, `-log(bhattacharyya_coefficient(...))`.
"""
bhattacharyya_distance(left, right; frequencies::Bool=true) =
    dissimilarity(Bhattacharyya(), left, right; frequencies)
bhattacharyya_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(Bhattacharyya(), data; frequencies, species)

"""
    kullback_leibler_divergence(left, right; frequencies=true, base=2, estimator=Plugin(), support=nothing)

Return Kullback-Leibler divergence ``D_{KL}(left \\Vert right)`` between
normalized abundance/probability vectors. This divergence is asymmetric and
returns `Inf` when `right` has zero probability where `left` has positive
probability.
"""
kullback_leibler_divergence(left, right; frequencies::Bool=true, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(KullbackLeibler(; base, estimator, support), left, right; frequencies)
kullback_leibler_divergence(data; frequencies::Bool=true, species=nothing, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(KullbackLeibler(; base, estimator, support), data; frequencies, species)

"""
    shannon_difference(left, right; frequencies=true, base=2, estimator=Plugin(), support=nothing)

Return the absolute difference between Shannon entropies of two assemblages.
"""
shannon_difference(left, right; frequencies::Bool=true, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(ShannonDifference(; base, estimator, support), left, right; frequencies)
shannon_difference(data; frequencies::Bool=true, species=nothing, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(ShannonDifference(; base, estimator, support), data; frequencies, species)

"""
    jensen_difference(left, right; frequencies=true, base=2, estimator=Plugin(), support=nothing)

Return the Jensen difference of Shannon entropy. This equals Jensen-Shannon
divergence for Shannon entropy.
"""
jensen_difference(left, right; frequencies::Bool=true, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(JensenDifference(; base, estimator, support), left, right; frequencies)
jensen_difference(data; frequencies::Bool=true, species=nothing, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(JensenDifference(; base, estimator, support), data; frequencies, species)

"""
    jensen_shannon_similarity(left, right; frequencies=true, base=2, distance=true, estimator=Plugin(), support=nothing)

Return one minus normalized Jensen-Shannon dissimilarity. By default the
normalization uses the square-root Jensen-Shannon distance; pass
`distance=false` to normalize the divergence instead.
"""
jensen_shannon_similarity(left, right; frequencies::Bool=true, base::Real=2,
        distance::Bool=true, estimator::ShannonEstimator=Plugin(), support=nothing) =
    similarity(JensenShannon(; base, distance, estimator, support), left, right; frequencies)
jensen_shannon_similarity(data; frequencies::Bool=true, species=nothing, base::Real=2,
        distance::Bool=true, estimator::ShannonEstimator=Plugin(), support=nothing) =
    similarity(JensenShannon(; base, distance, estimator, support), data; frequencies, species)

"""
    jensen_shannon_divergence(left, right; frequencies=true, base=2, estimator=Plugin(), support=nothing)

Return Jensen-Shannon divergence between probability vectors.
"""
jensen_shannon_divergence(left, right; frequencies::Bool=true, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(JensenShannon(; base, distance=false, estimator, support), left, right; frequencies)
jensen_shannon_divergence(data; frequencies::Bool=true, species=nothing, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(JensenShannon(; base, distance=false, estimator, support), data; frequencies, species)

"""
    jensen_shannon_distance(left, right; frequencies=true, base=2, estimator=Plugin(), support=nothing)

Return the square root of Jensen-Shannon divergence.
"""
jensen_shannon_distance(left, right; frequencies::Bool=true, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(JensenShannon(; base, distance=true, estimator, support), left, right; frequencies)
jensen_shannon_distance(data; frequencies::Bool=true, species=nothing, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing) =
    dissimilarity(JensenShannon(; base, distance=true, estimator, support), data; frequencies, species)

"""
    morisita_horn_similarity(left, right; frequencies=true)

Return Morisita-Horn abundance similarity.
"""
morisita_horn_similarity(left, right; frequencies::Bool=true) =
    similarity(MorisitaHorn(), left, right; frequencies)
morisita_horn_similarity(data; frequencies::Bool=true, species=nothing) =
    similarity(MorisitaHorn(), data; frequencies, species)

"""
    morisita_horn_distance(left, right; frequencies=true)

Return one minus Morisita-Horn similarity.
"""
morisita_horn_distance(left, right; frequencies::Bool=true) =
    dissimilarity(MorisitaHorn(), left, right; frequencies)
morisita_horn_distance(data; frequencies::Bool=true, species=nothing) =
    dissimilarity(MorisitaHorn(), data; frequencies, species)
