"""
Jaccard incidence similarity.

```math
J(A,B) = \\frac{|A \\cap B|}{|A \\cup B|}
```
"""
struct Jaccard <: DiversityIndex end

"""
Sorensen-Dice incidence similarity.

```math
S(A,B) = \\frac{2|A \\cap B|}{|A| + |B|}
```
"""
struct SorensenDice <: DiversityIndex end

"""
Overlap (Szymkiewicz-Simpson) incidence similarity.

```math
O(A,B) = \\frac{|A \\cap B|}{\\min(|A|, |B|)}
```
"""
struct Overlap <: DiversityIndex end

"""
Bray-Curtis abundance dissimilarity.

```math
BC(x,y) = \\frac{\\sum_i |x_i - y_i|}{\\sum_i (x_i + y_i)}
```
"""
struct BrayCurtis <: DiversityIndex end

"""
Ruzicka, or quantitative Jaccard, abundance similarity.

```math
R(x,y) = \\frac{\\sum_i \\min(x_i,y_i)}{\\sum_i \\max(x_i,y_i)}
```
"""
struct Ruzicka <: DiversityIndex end

"""
Total variation distance between probability vectors.

```math
TV(p,q) = \\frac{1}{2}\\sum_i |p_i - q_i|
```
"""
struct TotalVariation <: DiversityIndex end

"""
Manhattan, or L1, distance between probability vectors.

```math
d_1(p,q) = \\sum_i |p_i - q_i|
```
"""
struct Manhattan <: DiversityIndex end

"""
Euclidean, or L2, distance between probability vectors.

```math
d_2(p,q) = \\sqrt{\\sum_i (p_i - q_i)^2}
```
"""
struct Euclidean <: DiversityIndex end

"""
Averaged Canberra distance between abundance vectors.

```math
C(x,y) = \\frac{1}{m}\\sum_{i:x_i+y_i>0}\\frac{|x_i-y_i|}{x_i+y_i}
```
"""
struct Canberra <: DiversityIndex end

"""
Hellinger distance between probability vectors.

```math
H(p,q) = \\frac{1}{\\sqrt{2}}\\sqrt{\\sum_i (\\sqrt{p_i} - \\sqrt{q_i})^2}
```
"""
struct Hellinger <: DiversityIndex end

"""
Chord distance between square-root transformed probability vectors.

```math
d_c(p,q) = \\sqrt{\\sum_i (\\sqrt{p_i} - \\sqrt{q_i})^2}
```
"""
struct Chord <: DiversityIndex end

"""
Bhattacharyya coefficient and distance between probability vectors.

```math
BC(p,q) = \\sum_i \\sqrt{p_i q_i}, \\qquad
d_B(p,q) = -\\log BC(p,q)
```
"""
struct Bhattacharyya <: DiversityIndex end

"""
Jensen-Shannon divergence or distance between probability vectors.

`JensenShannon(; base=2, distance=true)` returns the square root of the
divergence from [`dissimilarity`](@ref). Set `distance=false` for the
divergence itself.
"""
struct JensenShannon <: DiversityIndex
    base::Float64
    distance::Bool
end

function JensenShannon(; base::Real=2, distance::Bool=true)
    base > 0 && base != 1 || throw(ArgumentError("base must be positive and not equal to 1"))
    return JensenShannon(float(base), distance)
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
struct MorisitaHorn <: DiversityIndex end

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

function _jensen_shannon_divergence(index::JensenShannon, left, right; frequencies::Bool=true)
    left_probability, right_probability = _aligned_probabilities(left, right; frequencies)
    mixture = (left_probability .+ right_probability) ./ 2
    return (_kl_divergence(left_probability, mixture, index.base) +
        _kl_divergence(right_probability, mixture, index.base)) / 2
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
    jensen_shannon_divergence(left, right; frequencies=true, base=2)

Return Jensen-Shannon divergence between probability vectors.
"""
jensen_shannon_divergence(left, right; frequencies::Bool=true, base::Real=2) =
    dissimilarity(JensenShannon(; base, distance=false), left, right; frequencies)
jensen_shannon_divergence(data; frequencies::Bool=true, species=nothing, base::Real=2) =
    dissimilarity(JensenShannon(; base, distance=false), data; frequencies, species)

"""
    jensen_shannon_distance(left, right; frequencies=true, base=2)

Return the square root of Jensen-Shannon divergence.
"""
jensen_shannon_distance(left, right; frequencies::Bool=true, base::Real=2) =
    dissimilarity(JensenShannon(; base, distance=true), left, right; frequencies)
jensen_shannon_distance(data; frequencies::Bool=true, species=nothing, base::Real=2) =
    dissimilarity(JensenShannon(; base, distance=true), data; frequencies, species)

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
