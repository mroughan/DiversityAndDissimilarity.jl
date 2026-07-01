"""Abstract supertype for Shannon entropy estimators."""
abstract type ShannonEstimator end

"""
Plugin estimator for Shannon entropy.

This is the maximum-likelihood estimator (MLE) obtained by substituting the
empirical probabilities into Shannon's definition.
"""
struct Plugin <: ShannonEstimator end

"""
Miller-Madow estimator for Shannon entropy.

For ``S`` observed positive-abundance species and sample size ``n``, the
correction added to Shannon entropy is

```math
\\frac{S - 1}{2n\\log b},
```

where ``b`` is the logarithm base.
"""
struct MillerMadow <: ShannonEstimator end

"""
Hausser-Strimmer shrinkage estimator for Shannon entropy.

This shrinks empirical probabilities towards the uniform distribution over the
known or observed support, then applies the plugin entropy estimator.
"""
struct HausserStrimmer <: ShannonEstimator end

"""
Basharin bias-reduced estimator for Shannon entropy.

For known support size ``S`` and sample size ``n``, the correction added to the
plugin estimator is

```math
\\frac{S - 1}{n\\log b},
```

where ``b`` is the logarithm base.
"""
struct Basharin <: ShannonEstimator end

"""
    AddGamma(gamma)

Add-gamma estimator for Shannon entropy.

With support size ``S`` and counts ``n_i``, probabilities are estimated as

```math
\\hat p_i = \\frac{n_i + \\gamma}{n + \\gamma S}.
```
"""
struct AddGamma <: ShannonEstimator
    gamma::Float64
    function AddGamma(gamma::Real)
        gamma >= 0 || throw(ArgumentError("gamma must be non-negative"))
        return new(float(gamma))
    end
end
AddGamma() = AddGamma(1)

"""
Chao-Shen estimator for Shannon entropy with possible unseen species.

This uses a Good-Turing coverage estimate and a Horvitz-Thompson adjustment.
It is intended for samples where the true support may be larger than the
observed support.
"""
struct ChaoShen <: ShannonEstimator end

"""Abstract supertype for diversity, similarity, and dissimilarity indices."""
abstract type DiversityIndex end

"""Abstract supertype for single-assemblage diversity, entropy, and richness indices."""
abstract type AlphaDiversityIndex <: DiversityIndex end

"""Species richness: the number of species/categories with positive abundance."""
struct Richness <: AlphaDiversityIndex end

"""
    Shannon(; base=2, estimator=Plugin())

Shannon entropy/index. Use `effective_diversity(Shannon(), x)` for the
corresponding Hill number, i.e. the effective number of species.

The default logarithm base is 2, so entropy is measured in bits unless another
base is provided. Use `estimator` to choose among [`Plugin`](@ref),
[`MillerMadow`](@ref), [`HausserStrimmer`](@ref), [`Basharin`](@ref),
[`AddGamma`](@ref), and [`ChaoShen`](@ref).

```math
H_b = -\\sum_i p_i \\log_b p_i
```
"""
struct Shannon{E<:ShannonEstimator} <: AlphaDiversityIndex
    base::Float64
    estimator::E
end
Shannon(; base::Real=2, estimator::ShannonEstimator=Plugin()) =
    Shannon(float(base), estimator)

"""
    Renyi(q; base=2)

Renyi entropy of order `q`.

The default logarithm base is 2, so entropy is measured in bits unless another
base is provided. At `q = 1`, Renyi entropy is evaluated as Shannon entropy.

```math
H_q = \\frac{1}{1-q}\\log_b\\left(\\sum_i p_i^q\\right)
```
"""
struct Renyi <: AlphaDiversityIndex
    q::Float64
    base::Float64
end
Renyi(q::Real; base::Real=2) = Renyi(float(q), float(base))

"""
    Tsallis(q; base=2)

Tsallis entropy of order `q`, scaled to the requested logarithm base.

The default logarithm base is 2. At `q = 1`, Tsallis entropy is evaluated as
Shannon entropy with the same base.

```math
T_q = \\frac{\\sum_i p_i^q - 1}{(1-q)\\log b}
```

!!! note "Non-standard base scaling"
    The ``\\log b`` denominator scales the result so that ``T_1`` equals Shannon
    entropy in the same base. This differs from the standard Tsallis definition
    ``(1 - \\sum_i p_i^q)/(q-1)``, which does not include a logarithm base factor.
    Values will differ from packages that use the standard definition by a factor
    of ``\\log b``.
"""
struct Tsallis <: AlphaDiversityIndex
    q::Float64
    base::Float64
end
Tsallis(q::Real; base::Real=2) = Tsallis(float(q), float(base))

"""
Simpson concentration.

```math
D = \\sum_i p_i^2
```
"""
struct Simpson <: AlphaDiversityIndex end

"""
Gini-Simpson diversity.

```math
1 - D = 1 - \\sum_i p_i^2
```
"""
struct GiniSimpson <: AlphaDiversityIndex end

"""
Greenberg's linguistic diversity index.

This is the probability that two randomly selected individuals have different
mother tongues. It is mathematically the same quantity as Gini-Simpson
diversity:

```math
1 - \\sum_i p_i^2
```

The index uses only relative mother-tongue frequencies. It does not account
for second-language use, language vitality, or linguistic distance between
languages.
"""
struct GreenbergDiversityIndex <: AlphaDiversityIndex end

"""
Linguistic diversity index (LDI).

This is an alias-by-convention for [`GreenbergDiversityIndex`](@ref), included
for linguistic-demography workflows. It returns the probability that two
randomly selected individuals have different mother tongues, equivalent to
[`GiniSimpson`](@ref).
"""
struct LinguisticDiversityIndex <: AlphaDiversityIndex end

"""
Inverse Simpson diversity.

```math
\\frac{1}{D} = \\frac{1}{\\sum_i p_i^2}
```
"""
struct InverseSimpson <: AlphaDiversityIndex end

"""
    Hill(q)

Hill diversity of order `q`. Orders 0, 1, and 2 correspond to richness,
Shannon effective diversity, and inverse Simpson diversity.

```math
{}^qD = \\left(\\sum_i p_i^q\\right)^{1/(1-q)}
```
"""
struct Hill <: AlphaDiversityIndex
    q::Float64
end
Hill(q::Real) = Hill(float(q))

"""
Chao1 asymptotic richness estimator.

For observed richness ``S_{obs}``, singleton count ``f_1``, and doubleton count
``f_2``, this package uses the bias-corrected abundance form

```math
\\hat S_{Chao1} = S_{obs} + \\frac{f_1(f_1 - 1)}{2(f_2 + 1)}.
```
"""
struct Chao1 <: AlphaDiversityIndex end

"""
    ACE(; threshold=10)

Abundance-based Coverage Estimator (ACE) for richness.

Species with counts up to `threshold` are treated as rare. The default
`threshold=10` follows the usual ACE convention.
"""
struct ACE <: AlphaDiversityIndex
    threshold::Int
    function ACE(threshold::Integer)
        threshold >= 1 || throw(ArgumentError("threshold must be positive"))
        return new(Int(threshold))
    end
end
ACE(; threshold::Integer=10) = ACE(threshold)

"""
Good-Turing sample coverage estimate.

For sample size ``n`` and singleton count ``f_1``, sample coverage is

```math
\\hat C = 1 - \\frac{f_1}{n}.
```
"""
struct SampleCoverage <: AlphaDiversityIndex end

"""
Pielou evenness, Shannon entropy divided by the maximum entropy for the
observed richness.

```math
J = \\frac{H}{\\log_b S}
```
"""
struct PielouEvenness <: AlphaDiversityIndex end

"""
Fisher's alpha diversity parameter.

Fisher's alpha solves

```math
S = \\alpha \\log\\left(1 + \\frac{n}{\\alpha}\\right),
```

where ``S`` is observed richness and ``n`` is total abundance.
"""
struct FisherAlpha <: AlphaDiversityIndex end

"""
    entropy(index, data; frequencies=true, support=nothing)

Evaluate an entropy index for `data`.

`entropy` is defined for [`Shannon`](@ref), [`Renyi`](@ref), and
[`Tsallis`](@ref). For dictionaries, values are abundances. Numeric vectors are
abundance vectors by default. Non-numeric vectors are raw observations. Pass
`frequencies=false` to treat a numeric vector as raw observations instead.

For [`Shannon`](@ref), `support` may be an integer support size or a collection
of known categories. Leave `support=nothing` when only the observed support is
known, or when using [`ChaoShen`](@ref) for possible unseen categories.

For community matrices, rows are samples and columns are taxa/categories. The
result is one entropy value per row.
"""
function entropy(index::Shannon, data; frequencies::Bool=true, support=nothing, species=nothing)
    if _is_table(data)
        return entropy(index, community_matrix(data; species); frequencies, support)
    end
    counts, n = _entropy_counts(data; frequencies, support)
    return _shannon_entropy(index.estimator, counts, n, index.base)
end

function entropy(index::Shannon{Plugin}, data::AbstractMatrix{<:Real}; frequencies::Bool=true, support=nothing, species=nothing)
    _check_matrix_frequencies(frequencies)
    if support === nothing
        _validate_community_matrix(data)
        return entropy(index, Validated(data); frequencies)
    end
    matrix_support = _matrix_support(support)
    return [entropy(index, row; frequencies, support=matrix_support) for row in _rows(data)]
end

function entropy(index::Shannon{Plugin}, v::Validated{<:AbstractMatrix{<:Real}}; frequencies::Bool=true, support=nothing, species=nothing)
    _check_matrix_frequencies(frequencies)
    if support !== nothing
        matrix_support = _matrix_support(support)
        return [entropy(index, row; frequencies, support=matrix_support) for row in _rows(v.data)]
    end
    data = v.data
    nsites = size(data, 1)
    row_totals = zeros(Float64, nsites)
    @inbounds for column in axes(data, 2)
        for row in axes(data, 1)
            row_totals[row] += data[row, column]
        end
    end
    result = zeros(Float64, nsites)
    log_base = log(index.base)
    @inbounds for column in axes(data, 2)
        for row in axes(data, 1)
            abundance = data[row, column]
            if abundance > 0
                probability = abundance / row_totals[row]
                result[row] -= probability * log(probability)
            end
        end
    end
    @inbounds for row in eachindex(result)
        result[row] /= log_base
    end
    return result
end

function entropy(index::Shannon, data::AbstractMatrix{<:Real}; frequencies::Bool=true, support=nothing, species=nothing)
    _check_matrix_frequencies(frequencies)
    matrix_support = _matrix_support(support)
    return [entropy(index, row; frequencies, support=matrix_support) for row in _rows(data)]
end

"""
    diversity(index, data; frequencies=true, support=nothing)

Evaluate a diversity index for `data`.

For dictionaries, values are abundances. Numeric vectors are abundance vectors
by default. Non-numeric vectors are raw observations. Pass `frequencies=false`
to treat a numeric vector as raw observations instead.

For entropy-family indices, `diversity` returns the corresponding effective
diversity. Use [`entropy`](@ref) when you want entropy units such as bits.
For community matrices, rows are samples and columns are taxa/categories. The
result is one diversity value per row.
"""
function diversity(::Richness, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return diversity(Richness(), community_matrix(data; species); frequencies)
    end
    return length(_abundances(data; frequencies))
end

function diversity(::Richness, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing)
    _check_matrix_frequencies(frequencies)
    _validate_community_matrix(data)
    return diversity(Richness(), Validated(data); frequencies)
end

function diversity(::Richness, v::Validated{<:AbstractMatrix{<:Real}}; frequencies::Bool=true, species=nothing)
    _check_matrix_frequencies(frequencies)
    data = v.data
    result = zeros(Int, size(data, 1))
    @inbounds for column in axes(data, 2)
        for row in axes(data, 1)
            result[row] += data[row, column] > 0
        end
    end
    return result
end

function entropy(index::Renyi, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return entropy(index, community_matrix(data; species); frequencies)
    end
    q = index.q
    if isapprox(q, 1.0; atol=sqrt(eps(Float64)))
        return entropy(Shannon(; base=index.base), data; frequencies)
    else
        p = proportions(data; frequencies)
        return log(sum(pi -> pi^q, p)) / ((1 - q) * log(index.base))
    end
end

entropy(index::Renyi, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [entropy(index, row; frequencies) for row in _rows(data)])

function entropy(index::Tsallis, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return entropy(index, community_matrix(data; species); frequencies)
    end
    q = index.q
    if isapprox(q, 1.0; atol=sqrt(eps(Float64)))
        return entropy(Shannon(; base=index.base), data; frequencies)
    else
        p = proportions(data; frequencies)
        return (sum(pi -> pi^q, p) - 1) / ((1 - q) * log(index.base))
    end
end

entropy(index::Tsallis, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [entropy(index, row; frequencies) for row in _rows(data)])

function diversity(index::Shannon, data; frequencies::Bool=true, support=nothing, species=nothing)
    if _is_table(data)
        return diversity(index, community_matrix(data; species); frequencies, support)
    end
    return effective_diversity(index, data; frequencies, support)
end

function diversity(index::Shannon, data::AbstractMatrix{<:Real}; frequencies::Bool=true, support=nothing, species=nothing)
    _check_matrix_frequencies(frequencies)
    return effective_diversity(index, data; frequencies, support)
end

diversity(index::Renyi, data; frequencies::Bool=true, species=nothing) =
    effective_diversity(index, _community_input(data; species); frequencies)
diversity(index::Renyi, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); effective_diversity(index, data; frequencies))

diversity(index::Tsallis, data; frequencies::Bool=true, species=nothing) =
    effective_diversity(index, _community_input(data; species); frequencies)
diversity(index::Tsallis, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); effective_diversity(index, data; frequencies))

diversity(::Simpson, data; frequencies::Bool=true, species=nothing) =
    _is_table(data) ?
        diversity(Simpson(), community_matrix(data; species); frequencies) :
        sum(abs2, proportions(data; frequencies, species))
diversity(index::Simpson, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

diversity(::GiniSimpson, data; frequencies::Bool=true, species=nothing) =
    _is_table(data) ?
        diversity(GiniSimpson(), community_matrix(data; species); frequencies) :
        1 - diversity(Simpson(), data; frequencies, species)
diversity(index::GiniSimpson, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

diversity(::Union{GreenbergDiversityIndex,LinguisticDiversityIndex}, data; frequencies::Bool=true, species=nothing) =
    _is_table(data) ?
        diversity(GiniSimpson(), community_matrix(data; species); frequencies) :
        diversity(GiniSimpson(), data; frequencies, species)
diversity(index::Union{GreenbergDiversityIndex,LinguisticDiversityIndex}, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

diversity(::InverseSimpson, data; frequencies::Bool=true, species=nothing) =
    _is_table(data) ?
        diversity(InverseSimpson(), community_matrix(data; species); frequencies) :
        inv(diversity(Simpson(), data; frequencies, species))
diversity(index::InverseSimpson, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

function diversity(index::Hill, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return diversity(index, community_matrix(data; species); frequencies)
    end
    q = index.q
    p = proportions(data; frequencies)
    if isapprox(q, 0.0; atol=eps(Float64))
        return length(p)
    elseif isapprox(q, 1.0; atol=sqrt(eps(Float64)))
        return diversity(Shannon(), data; frequencies)
    else
        return sum(pi -> pi^q, p)^(1 / (1 - q))
    end
end

diversity(index::Hill, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

function diversity(::Chao1, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return diversity(Chao1(), community_matrix(data; species); frequencies)
    end
    abundance = _abundances(data; frequencies)
    f1 = count(==(1), abundance)
    f2 = count(==(2), abundance)
    return length(abundance) + float(f1) * (f1 - 1) / (2 * (f2 + 1))
end

diversity(index::Chao1, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

function diversity(index::ACE, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return diversity(index, community_matrix(data; species); frequencies)
    end
    abundance = _abundances(data; frequencies)
    threshold = index.threshold
    rare = [value for value in abundance if value <= threshold]
    abundant = count(>(threshold), abundance)
    isempty(rare) && return float(abundant)

    n_rare = sum(rare)
    f1 = count(==(1), rare)
    coverage = 1 - f1 / n_rare
    coverage > 0 || throw(ArgumentError(
        "ACE coverage estimate is zero: all $f1 rare-species counts are singletons " *
        "(f₁=$f1, n_rare=$n_rare); try a larger dataset or a lower threshold"))

    s_rare = length(rare)
    frequency_counts = [count(==(i), rare) for i in 1:threshold]
    gamma_squared = if n_rare <= 1
        0.0
    else
        numerator = sum(i -> i * (i - 1) * frequency_counts[i], 1:threshold)
        max(s_rare / coverage * numerator / (n_rare * (n_rare - 1)) - 1, 0.0)
    end
    return abundant + s_rare / coverage + f1 / coverage * gamma_squared
end

diversity(index::ACE, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

function diversity(::SampleCoverage, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return diversity(SampleCoverage(), community_matrix(data; species); frequencies)
    end
    abundance = _abundances(data; frequencies)
    n = sum(abundance)
    f1 = count(==(1), abundance)
    return 1 - f1 / n
end

diversity(index::SampleCoverage, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

function diversity(::PielouEvenness, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return diversity(PielouEvenness(), community_matrix(data; species); frequencies)
    end
    observed_richness = richness(data; frequencies)
    observed_richness > 1 || return 1.0
    return entropy(Shannon(; base=ℯ), data; frequencies) / log(observed_richness)
end

diversity(index::PielouEvenness, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

function diversity(::FisherAlpha, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return diversity(FisherAlpha(), community_matrix(data; species); frequencies)
    end
    abundance = _abundances(data; frequencies)
    observed_richness = length(abundance)
    total = sum(abundance)
    observed_richness <= total || throw(ArgumentError("observed richness cannot exceed total abundance"))
    observed_richness == total && return Inf
    observed_richness == 1 && total == 1 && return Inf
    return _fisher_alpha(observed_richness, total)
end

diversity(index::FisherAlpha, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [diversity(index, row; frequencies) for row in _rows(data)])

function _fisher_alpha(observed_richness, total)
    target = float(observed_richness)
    n = float(total)
    lower = 0.0
    upper = max(target, 1.0)
    while upper * log1p(n / upper) < target
        upper *= 2
    end
    for _ in 1:100
        middle = (lower + upper) / 2
        value = middle * log1p(n / middle)
        if value < target
            lower = middle
        else
            upper = middle
        end
    end
    return (lower + upper) / 2
end

"""
    effective_diversity(index, data; frequencies=true)

Return an index as an effective number of species when the transformation is
standard. For Shannon entropy this is

```math
{}^1D = b^H.
```

For community matrices, rows are samples and columns are taxa/categories. The
result is one effective-diversity value per row.
"""
effective_diversity(index::Hill, data; frequencies::Bool=true, species=nothing) =
    diversity(index, data; frequencies, species)
effective_diversity(::Richness, data; frequencies::Bool=true, species=nothing) =
    diversity(Richness(), data; frequencies, species)
effective_diversity(index::Shannon, data; frequencies::Bool=true, support=nothing, species=nothing) =
    index.base ^ entropy(index, data; frequencies, support, species)
function effective_diversity(index::Shannon, data::AbstractMatrix{<:Real}; frequencies::Bool=true, support=nothing, species=nothing)
    _check_matrix_frequencies(frequencies)
    matrix_support = _matrix_support(support)
    return [effective_diversity(index, row; frequencies, support=matrix_support) for row in _rows(data)]
end
effective_diversity(index::Renyi, data; frequencies::Bool=true, species=nothing) =
    index.base ^ entropy(index, data; frequencies, species)
effective_diversity(index::Renyi, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [effective_diversity(index, row; frequencies) for row in _rows(data)])
function effective_diversity(index::Tsallis, data; frequencies::Bool=true, species=nothing)
    if _is_table(data)
        return effective_diversity(index, community_matrix(data; species); frequencies)
    end
    q = index.q
    if isapprox(q, 1.0; atol=sqrt(eps(Float64)))
        return effective_diversity(Shannon(; base=index.base), data; frequencies)
    else
        value = entropy(index, data; frequencies)
        return (1 + (1 - q) * log(index.base) * value)^(1 / (1 - q))
    end
end
effective_diversity(index::Tsallis, data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing) =
    (_check_matrix_frequencies(frequencies); [effective_diversity(index, row; frequencies) for row in _rows(data)])
effective_diversity(::Simpson, data; frequencies::Bool=true, species=nothing) =
    diversity(InverseSimpson(), data; frequencies, species)
effective_diversity(::InverseSimpson, data; frequencies::Bool=true, species=nothing) =
    diversity(InverseSimpson(), data; frequencies, species)

_matrix_support(support) = support === nothing || support isa Integer ? support : length(collect(support))

function _entropy_counts(data; frequencies::Bool=true, support=nothing)
    observed = counts(data)
    if data isa AbstractVector{<:Real} && frequencies
        for value in data
            value < 0 && throw(ArgumentError("abundances must be non-negative"))
            isfinite(value) || throw(ArgumentError("abundances must be finite"))
        end
        observed = Dict(i => float(value) for (i, value) in pairs(data))
    end

    categories = collect(keys(observed))
    if support === nothing
        abundance_values = collect(float(v) for v in Base.values(observed))
    elseif support isa Integer
        support >= length(categories) || throw(ArgumentError("support size must be at least the observed support size"))
        abundance_values = collect(float(v) for v in Base.values(observed))
        append!(abundance_values, zeros(Float64, support - length(categories)))
    else
        known = collect(support)
        observed_set = Set(categories)
        known_set = Set(known)
        observed_set ⊆ known_set || throw(ArgumentError("support must include every observed category"))
        abundance_values = [float(get(observed, category, 0)) for category in known]
    end

    any(v -> v < 0 || !isfinite(v), abundance_values) && throw(ArgumentError("abundances must be non-negative and finite"))
    n = sum(abundance_values)
    n > 0 || throw(ArgumentError("total abundance must be positive"))
    return (abundance_values, n)
end

_positive_probabilities(counts, n) = [count / n for count in counts if count > 0]

function _plugin_entropy_from_probabilities(probabilities, base)
    return -sum(p -> p > 0 ? p * log(p) : 0.0, probabilities) / log(base)
end

function _shannon_entropy(::Plugin, counts, n, base)
    return _plugin_entropy_from_probabilities(_positive_probabilities(counts, n), base)
end

function _shannon_entropy(::MillerMadow, counts, n, base)
    m_positive = count(>(0), counts)
    return _shannon_entropy(Plugin(), counts, n, base) + (m_positive - 1) / (2n * log(base))
end

function _shannon_entropy(::Basharin, counts, n, base)
    support_size = length(counts)
    return _shannon_entropy(Plugin(), counts, n, base) + (support_size - 1) / (n * log(base))
end

function _shannon_entropy(estimator::AddGamma, counts, n, base)
    support_size = length(counts)
    denominator = n + estimator.gamma * support_size
    probabilities = [(count + estimator.gamma) / denominator for count in counts]
    return _plugin_entropy_from_probabilities(probabilities, base)
end

function _shannon_entropy(::HausserStrimmer, counts, n, base)
    support_size = length(counts)
    target = inv(support_size)
    probabilities = [count / n for count in counts]
    denominator = sum(p -> (target - p)^2, probabilities)
    shrinkage = if n <= 1 || denominator <= eps(Float64)
        1.0
    else
        numerator = sum(p -> p * (1 - p), probabilities) / (n - 1)
        clamp(numerator / denominator, 0.0, 1.0)
    end
    shrinkage_probabilities = [shrinkage * target + (1 - shrinkage) * p for p in probabilities]
    return _plugin_entropy_from_probabilities(shrinkage_probabilities, base)
end

function _shannon_entropy(::ChaoShen, counts, n, base)
    positive_counts = [count for count in counts if count > 0]
    singletons = count(==(1), positive_counts)
    coverage = 1 - singletons / n
    coverage > 0 || throw(ArgumentError(
        "Chao-Shen coverage estimate is zero: all $n observations are singletons " *
        "(singletons=$singletons, n=$n); try a larger dataset or Plugin/MillerMadow estimator"))
    probabilities = [coverage * count / n for count in positive_counts]
    terms = map(probabilities) do p
        detection_probability = 1 - (1 - p)^n
        detection_probability > 0 ? p * log(p) / detection_probability : 0.0
    end
    return -sum(terms) / log(base)
end

"""
    entropy_variance(index, data; frequencies=true, support=nothing)

Estimate the variance of a Shannon entropy estimator.

Analytic variance estimates are available for [`Plugin`](@ref),
[`MillerMadow`](@ref), [`Basharin`](@ref), and [`ChaoShen`](@ref). For
[`ChaoShen`](@ref), this is an approximate Horvitz-Thompson detection variance.
Use [`bootstrap`](@ref) or [`jackknife`](@ref) for estimators without an analytic
variance estimate.
"""
function entropy_variance(index::Shannon, data; frequencies::Bool=true, support=nothing, species=nothing)
    if _is_table(data)
        return entropy_variance(index, community_matrix(data; species); frequencies, support)
    end
    counts, n = _entropy_counts(data; frequencies, support)
    return _shannon_variance(index.estimator, counts, n, index.base)
end

function entropy_variance(index::Shannon, data::AbstractMatrix{<:Real}; frequencies::Bool=true, support=nothing, species=nothing)
    _check_matrix_frequencies(frequencies)
    matrix_support = _matrix_support(support)
    return [entropy_variance(index, row; frequencies, support=matrix_support) for row in _rows(data)]
end

_shannon_variance(::Plugin, counts, n, base) =
    _plugin_shannon_variance(counts, n, base; second_order=false)
_shannon_variance(::MillerMadow, counts, n, base) =
    _plugin_shannon_variance(counts, n, base; second_order=false)
_shannon_variance(::Basharin, counts, n, base) =
    _plugin_shannon_variance(counts, n, base; second_order=true)

function _shannon_variance(estimator::Union{HausserStrimmer,AddGamma}, counts, n, base)
    throw(ArgumentError("analytic variance is not implemented for $(typeof(estimator)); use bootstrap or jackknife"))
end

function _plugin_shannon_variance(counts, n, base; second_order::Bool=false)
    p = _positive_probabilities(counts, n)
    logp = log.(p)
    first_order = (sum(pi * li^2 for (pi, li) in zip(p, logp)) - sum(pi * li for (pi, li) in zip(p, logp))^2) / n
    correction = second_order ? (count(>(0), counts) - 1) / (2n^2) : 0.0
    return max((first_order + correction) / log(base)^2, 0.0)
end

function _shannon_variance(::ChaoShen, counts, n, base)
    positive_counts = [count for count in counts if count > 0]
    singletons = count(==(1), positive_counts)
    coverage = 1 - singletons / n
    coverage > 0 || throw(ArgumentError(
        "Chao-Shen coverage estimate is zero: all $n observations are singletons " *
        "(singletons=$singletons, n=$n); try a larger dataset or Plugin/MillerMadow estimator"))
    variance = 0.0
    for count_value in positive_counts
        probability = coverage * count_value / n
        detection_probability = 1 - (1 - probability)^n
        if detection_probability > 0
            contribution = probability * log(probability)
            variance += (1 - detection_probability) * contribution^2 / detection_probability^2
        end
    end
    return max(variance / log(base)^2, 0.0)
end

"""
    entropy_confint(index, data; level=0.95, frequencies=true, support=nothing)

Return a normal-approximation confidence interval for Shannon entropy using
[`entropy_variance`](@ref).
"""
function entropy_confint(index::Shannon, data; level::Real=0.95, frequencies::Bool=true, support=nothing, species=nothing)
    if _is_table(data)
        return entropy_confint(index, community_matrix(data; species); level, frequencies, support)
    end
    estimate = entropy(index, data; frequencies, support)
    variance = entropy_variance(index, data; frequencies, support)
    return _normal_interval(estimate, variance, level)
end

function entropy_confint(index::Shannon, data::AbstractMatrix{<:Real}; level::Real=0.95, frequencies::Bool=true, support=nothing, species=nothing)
    _check_matrix_frequencies(frequencies)
    matrix_support = _matrix_support(support)
    return [entropy_confint(index, row; level, frequencies, support=matrix_support) for row in _rows(data)]
end

function _normal_interval(estimate, variance, level)
    _check_confidence_level(level)
    stderr = sqrt(max(variance, 0.0))
    z = _normal_quantile((1 + float(level)) / 2)
    return (
        estimate=estimate,
        variance=variance,
        stderr=stderr,
        level=float(level),
        lower=estimate - z * stderr,
        upper=estimate + z * stderr,
    )
end

function _check_confidence_level(level)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    return nothing
end

# Acklam's rational approximation for the standard-normal quantile.
function _normal_quantile(p::Real)
    0 < p < 1 || throw(ArgumentError("probability must be between 0 and 1"))
    a = (-3.969683028665376e1, 2.209460984245205e2, -2.759285104469687e2,
        1.383577518672690e2, -3.066479806614716e1, 2.506628277459239)
    b = (-5.447609879822406e1, 1.615858368580409e2, -1.556989798598866e2,
        6.680131188771972e1, -1.328068155288572e1)
    c = (-7.784894002430293e-3, -3.223964580411365e-1, -2.400758277161838,
        -2.549732539343734, 4.374664141464968, 2.938163982698783)
    d = (7.784695709041462e-3, 3.224671290700398e-1, 2.445134137142996,
        3.754408661907416)
    plow = 0.02425
    phigh = 1 - plow
    if p < plow
        q = sqrt(-2log(p))
        return (((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
            ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1)
    elseif p <= phigh
        q = p - 0.5
        r = q^2
        return (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
            (((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1)
    else
        q = sqrt(-2log(1 - p))
        return -(((((c[1] * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) * q + c[6]) /
            ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1)
    end
end

"""
    bootstrap(index, data; nboot=1000, level=0.95, quantity=:entropy, rng=nothing)

Bootstrap a Shannon estimator from empirical abundance counts.

`quantity=:entropy` bootstraps entropy values. `quantity=:diversity` bootstraps
effective Shannon diversity values. The returned named tuple contains the
estimate, bootstrap variance, standard error, percentile interval, and
replicates.
"""
function bootstrap(index::Shannon, data; nboot::Integer=1000, level::Real=0.95, quantity::Symbol=:entropy,
        frequencies::Bool=true, support=nothing, species=nothing, rng=nothing)
    if _is_table(data)
        return bootstrap(index, community_matrix(data; species); nboot, level, quantity, frequencies, support, rng)
    end
    nboot > 1 || throw(ArgumentError("nboot must be greater than 1"))
    _check_confidence_level(level)
    counts, n = _integer_entropy_counts(data; frequencies, support)
    estimate = _shannon_quantity(index, counts, n, quantity)
    replicates = Vector{Float64}(undef, nboot)
    for replicate in 1:nboot
        boot_counts = _multinomial_resample_counts(rng, counts, n)
        replicates[replicate] = _shannon_quantity(index, boot_counts, n, quantity)
    end
    return _resampling_summary(estimate, replicates, level)
end

function bootstrap(index::Shannon, data::AbstractMatrix{<:Real}; nboot::Integer=1000, level::Real=0.95,
        quantity::Symbol=:entropy, frequencies::Bool=true, support=nothing, species=nothing,
        rng=nothing)
    _check_matrix_frequencies(frequencies)
    matrix_support = _matrix_support(support)
    return [bootstrap(index, row; nboot, level, quantity, frequencies, support=matrix_support, rng) for row in _rows(data)]
end

"""
    alpha_diversity(data; frequencies=true, species=nothing, base=2, estimator=Plugin(),
                    support=nothing, threshold=10)

Return a compact alpha-diversity summary for common exploratory workflows.

For one assemblage, the result is a named tuple containing observed richness,
Shannon entropy, Shannon effective diversity, Simpson concentration,
Gini-Simpson diversity, inverse Simpson diversity, Chao1, ACE, and sample
coverage. For community matrices and Tables.jl-compatible inputs, the result is
one named tuple per row/sample.
"""
function alpha_diversity(data; frequencies::Bool=true, species=nothing, base::Real=2,
        estimator::ShannonEstimator=Plugin(), support=nothing, threshold::Integer=10)
    if _is_table(data)
        return alpha_diversity(
            community_matrix(data; species);
            frequencies,
            base,
            estimator,
            support,
            threshold,
        )
    end
    return (
        richness=richness(data; frequencies),
        shannon_entropy=shannon_entropy(data; base, estimator, frequencies, support),
        shannon_diversity=shannon_diversity(data; base, estimator, frequencies, support),
        simpson=simpson_index(data; frequencies),
        gini_simpson=gini_simpson_index(data; frequencies),
        inverse_simpson=inverse_simpson_index(data; frequencies),
        chao1=chao1(data; frequencies),
        ace=ace(data; frequencies, threshold),
        sample_coverage=sample_coverage(data; frequencies),
        pielou_evenness=pielou_evenness(data; frequencies),
        fisher_alpha=fisher_alpha(data; frequencies),
    )
end

function alpha_diversity(data::AbstractMatrix{<:Real}; frequencies::Bool=true, species=nothing,
        base::Real=2, estimator::ShannonEstimator=Plugin(), support=nothing,
        threshold::Integer=10)
    _check_matrix_frequencies(frequencies)
    if estimator isa Plugin && support === nothing
        _validate_community_matrix(data)
        return _alpha_diversity_matrix_plugin(data, float(base), Int(threshold))
    end
    matrix_support = _matrix_support(support)
    return [
        alpha_diversity(
            row;
            frequencies,
            base,
            estimator,
            support=matrix_support,
            threshold,
        )
        for row in _rows(data)
    ]
end

function alpha_diversity(v::Validated{<:AbstractMatrix{<:Real}}; frequencies::Bool=true, species=nothing,
        base::Real=2, estimator::ShannonEstimator=Plugin(), support=nothing,
        threshold::Integer=10)
    _check_matrix_frequencies(frequencies)
    if estimator isa Plugin && support === nothing
        return _alpha_diversity_matrix_plugin(v.data, float(base), Int(threshold))
    end
    matrix_support = _matrix_support(support)
    return [
        alpha_diversity(row; frequencies, base, estimator, support=matrix_support, threshold)
        for row in _rows(v.data)
    ]
end

function _alpha_diversity_matrix_plugin(data::AbstractMatrix{<:Real}, base::Float64, threshold::Int)
    threshold >= 1 || throw(ArgumentError("threshold must be positive"))
    result = Vector{NamedTuple{
        (:richness, :shannon_entropy, :shannon_diversity, :simpson,
            :gini_simpson, :inverse_simpson, :chao1, :ace, :sample_coverage,
            :pielou_evenness, :fisher_alpha),
        Tuple{Int,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64},
    }}(undef, size(data, 1))
    log_base = log(base)
    nsites = size(data, 1)
    row_totals = zeros(Float64, nsites)
    richness = zeros(Int, nsites)
    f1 = zeros(Int, nsites)
    f2 = zeros(Int, nsites)
    rare_count = zeros(Int, nsites)
    abundant = zeros(Int, nsites)
    n_rare = zeros(Float64, nsites)
    ace_numerator = zeros(Float64, nsites)
    shannon_values = zeros(Float64, nsites)
    simpson_values = zeros(Float64, nsites)

    @inbounds for column in axes(data, 2)
        for row in axes(data, 1)
            abundance = float(data[row, column])
            row_totals[row] += abundance
            if abundance > 0
                richness[row] += 1
                f1[row] += abundance == 1
                f2[row] += abundance == 2
                if abundance <= threshold
                    rare_count[row] += 1
                    n_rare[row] += abundance
                    if abundance == floor(abundance)
                        ace_numerator[row] += abundance * (abundance - 1)
                    end
                else
                    abundant[row] += 1
                end
            end
        end
    end

    @inbounds for column in axes(data, 2)
        for row in axes(data, 1)
            abundance = float(data[row, column])
            if abundance > 0
                probability = abundance / row_totals[row]
                shannon_values[row] -= probability * log(probability)
                simpson_values[row] += probability^2
            end
        end
    end

    @inbounds for row in axes(data, 1)
        shannon_value = shannon_values[row] / log_base
        chao1_value = richness[row] + float(f1[row]) * (f1[row] - 1) / (2 * (f2[row] + 1))
        ace_value = if rare_count[row] == 0
            float(abundant[row])
        else
            coverage = 1 - f1[row] / n_rare[row]
            coverage > 0 || throw(ArgumentError(
                "ACE coverage estimate is zero at row $row: all $(f1[row]) rare-species counts are singletons " *
                "(f₁=$(f1[row]), n_rare=$(n_rare[row])); try a larger dataset or a lower threshold"))
            gamma_squared = if n_rare[row] <= 1
                0.0
            else
                max(rare_count[row] / coverage * ace_numerator[row] / (n_rare[row] * (n_rare[row] - 1)) - 1, 0.0)
            end
            abundant[row] + rare_count[row] / coverage + f1[row] / coverage * gamma_squared
        end

        richness[row] <= row_totals[row] || throw(ArgumentError("observed richness cannot exceed total abundance"))
        fisher_value = richness[row] == row_totals[row] ? Inf : _fisher_alpha(richness[row], row_totals[row])
        result[row] = (
            richness=richness[row],
            shannon_entropy=shannon_value,
            shannon_diversity=base^shannon_value,
            simpson=simpson_values[row],
            gini_simpson=1 - simpson_values[row],
            inverse_simpson=inv(simpson_values[row]),
            chao1=chao1_value,
            ace=ace_value,
            sample_coverage=1 - f1[row] / row_totals[row],
            pielou_evenness=richness[row] > 1 ? (shannon_values[row] / log(richness[row])) : 1.0,
            fisher_alpha=fisher_value,
        )
    end
    return result
end

"""
    jackknife(index, data; level=0.95, quantity=:entropy)

Delete-one jackknife for Shannon entropy or Shannon effective diversity.
"""
function jackknife(index::Shannon, data; level::Real=0.95, quantity::Symbol=:entropy,
        frequencies::Bool=true, support=nothing, species=nothing)
    if _is_table(data)
        return jackknife(index, community_matrix(data; species); level, quantity, frequencies, support)
    end
    _check_confidence_level(level)
    counts, n = _integer_entropy_counts(data; frequencies, support)
    n > 1 || throw(ArgumentError("jackknife requires sample size greater than one"))
    estimate = _shannon_quantity(index, counts, n, quantity)
    leave_values = Float64[]
    weights = Float64[]
    for (i, count_value) in pairs(counts)
        count_value <= 0 && continue
        reduced = copy(counts)
        reduced[i] -= 1
        push!(leave_values, _shannon_quantity(index, reduced, n - 1, quantity))
        push!(weights, count_value)
    end
    mean_leave = sum(w * value for (w, value) in zip(weights, leave_values)) / n
    variance = (n - 1) / n * sum(w * (value - mean_leave)^2 for (w, value) in zip(weights, leave_values))
    bias_corrected = n * estimate - (n - 1) * mean_leave
    interval = _normal_interval(bias_corrected, variance, level)
    return (
        estimate=estimate,
        bias_corrected=bias_corrected,
        variance=variance,
        stderr=interval.stderr,
        level=interval.level,
        lower=interval.lower,
        upper=interval.upper,
        leave_one_out=leave_values,
    )
end

function jackknife(index::Shannon, data::AbstractMatrix{<:Real}; level::Real=0.95, quantity::Symbol=:entropy,
        frequencies::Bool=true, support=nothing, species=nothing)
    _check_matrix_frequencies(frequencies)
    matrix_support = _matrix_support(support)
    return [jackknife(index, row; level, quantity, frequencies, support=matrix_support) for row in _rows(data)]
end

function _integer_entropy_counts(data; frequencies::Bool=true, support=nothing)
    raw_counts, _ = _entropy_counts(data; frequencies, support)
    integer_counts = Int[]
    for value in raw_counts
        rounded = round(Int, value)
        isapprox(value, rounded; atol=sqrt(eps(Float64))) || throw(ArgumentError("bootstrap and jackknife require integer count data"))
        push!(integer_counts, rounded)
    end
    n = sum(integer_counts)
    n > 0 || throw(ArgumentError("total abundance must be positive"))
    return integer_counts, n
end

function _multinomial_resample_counts(rng, counts, n)
    probabilities = [count / n for count in counts]
    cumulative = cumsum(probabilities)
    cumulative[end] = 1.0
    result = zeros(Int, length(counts))
    for _ in 1:n
        draw = rng === nothing ? rand() : rand(rng)
        index = searchsortedfirst(cumulative, draw)
        result[index] += 1
    end
    return result
end

function _shannon_quantity(index, counts, n, quantity)
    entropy_value = _shannon_entropy(index.estimator, counts, n, index.base)
    if quantity === :entropy
        return entropy_value
    elseif quantity === :diversity
        return index.base^entropy_value
    else
        throw(ArgumentError("quantity must be :entropy or :diversity"))
    end
end

function _resampling_summary(estimate, replicates, level)
    variance = _sample_variance(replicates)
    alpha = (1 - float(level)) / 2
    return (
        estimate=estimate,
        variance=variance,
        stderr=sqrt(max(variance, 0.0)),
        level=float(level),
        lower=_sample_quantile(replicates, alpha),
        upper=_sample_quantile(replicates, 1 - alpha),
        replicates=replicates,
    )
end

function _sample_variance(values)
    n = length(values)
    n > 1 || throw(ArgumentError("at least two values are required"))
    mean_value = sum(values) / n
    return sum(value -> (value - mean_value)^2, values) / (n - 1)
end

function _sample_quantile(values, probability)
    sorted = sort(collect(values))
    n = length(sorted)
    n > 0 || throw(ArgumentError("at least one value is required"))
    position = 1 + (n - 1) * probability
    lower_index = floor(Int, position)
    upper_index = ceil(Int, position)
    lower_index == upper_index && return sorted[lower_index]
    weight = position - lower_index
    return (1 - weight) * sorted[lower_index] + weight * sorted[upper_index]
end

"""
    richness(data; frequencies=true)

Return species richness for `data`.
"""
richness(data; frequencies::Bool=true, species=nothing) =
    diversity(Richness(), data; frequencies, species)

"""
    chao1(data; frequencies=true)

Return the Chao1 richness estimate for `data`.
"""
chao1(data; frequencies::Bool=true, species=nothing) =
    diversity(Chao1(), data; frequencies, species)

"""
    ace(data; frequencies=true, threshold=10)

Return the abundance-based coverage estimator (ACE) for richness.
"""
ace(data; frequencies::Bool=true, threshold::Integer=10, species=nothing) =
    diversity(ACE(; threshold), data; frequencies, species)

"""
    sample_coverage(data; frequencies=true)

Return the Good-Turing sample coverage estimate for `data`.
"""
sample_coverage(data; frequencies::Bool=true, species=nothing) =
    diversity(SampleCoverage(), data; frequencies, species)

"""
    pielou_evenness(data; frequencies=true)

Return Pielou evenness, Shannon entropy divided by the maximum entropy for the
observed richness.
"""
pielou_evenness(data; frequencies::Bool=true, species=nothing) =
    diversity(PielouEvenness(), data; frequencies, species)

"""
    fisher_alpha(data; frequencies=true)

Return Fisher's alpha diversity parameter.
"""
fisher_alpha(data; frequencies::Bool=true, species=nothing) =
    diversity(FisherAlpha(), data; frequencies, species)

"""
    shannon_entropy(data; base=2, estimator=Plugin(), frequencies=true, support=nothing)

Return Shannon entropy for `data`.

The default logarithm base is 2, so entropy is measured in bits. Pass `support`
when using finite-support estimators and the support is known.
"""
shannon_entropy(data; base::Real=2, estimator::ShannonEstimator=Plugin(), frequencies::Bool=true, support=nothing, species=nothing) =
    entropy(Shannon(; base, estimator), data; frequencies, support, species)

"""
    shannon_variance(data; base=2, estimator=Plugin(), frequencies=true, support=nothing)

Return the estimated variance of a Shannon entropy estimator.
"""
shannon_variance(data; base::Real=2, estimator::ShannonEstimator=Plugin(), frequencies::Bool=true, support=nothing, species=nothing) =
    entropy_variance(Shannon(; base, estimator), data; frequencies, support, species)

"""
    shannon_confint(data; base=2, estimator=Plugin(), level=0.95, frequencies=true, support=nothing)

Return a normal-approximation confidence interval for Shannon entropy.
"""
shannon_confint(data; base::Real=2, estimator::ShannonEstimator=Plugin(), level::Real=0.95,
        frequencies::Bool=true, support=nothing, species=nothing) =
    entropy_confint(Shannon(; base, estimator), data; level, frequencies, support, species)

"""
    shannon(data; base=2, estimator=Plugin(), frequencies=true, support=nothing)

Alias for [`shannon_entropy`](@ref).
"""
shannon(data; base::Real=2, estimator::ShannonEstimator=Plugin(), frequencies::Bool=true, support=nothing, species=nothing) =
    shannon_entropy(data; base, estimator, frequencies, support, species)

"""
    shannon_diversity(data; base=2, estimator=Plugin(), frequencies=true, support=nothing)

Return Shannon effective diversity for `data`.

```math
{}^1D = b^H
```

The default logarithm base is 2, matching [`Shannon`](@ref).
"""
shannon_diversity(data; base::Real=2, estimator::ShannonEstimator=Plugin(), frequencies::Bool=true, support=nothing, species=nothing) =
    effective_diversity(Shannon(; base, estimator), data; frequencies, support, species)

"""
    renyi_entropy(data, q; base=2, frequencies=true)

Return Renyi entropy of order `q` for `data`.

The default logarithm base is 2, matching [`Renyi`](@ref).

```math
H_q = \\frac{1}{1-q}\\log_b\\left(\\sum_i p_i^q\\right)
```
"""
renyi_entropy(data, q::Real; base::Real=2, frequencies::Bool=true, species=nothing) =
    entropy(Renyi(q; base), data; frequencies, species)

"""
    renyi(data, q; base=2, frequencies=true)

Alias for [`renyi_entropy`](@ref).
"""
renyi(data, q::Real; base::Real=2, frequencies::Bool=true, species=nothing) =
    renyi_entropy(data, q; base, frequencies, species)

"""
    renyi_diversity(data, q; base=2, frequencies=true)

Return the effective diversity corresponding to Renyi entropy of order `q`.

```math
{}^qD = b^{H_q}
```
"""
renyi_diversity(data, q::Real; base::Real=2, frequencies::Bool=true, species=nothing) =
    effective_diversity(Renyi(q; base), data; frequencies, species)

"""
    tsallis_entropy(data, q; base=2, frequencies=true)

Return Tsallis entropy of order `q` for `data`.

The default logarithm base is 2, matching [`Tsallis`](@ref).

```math
T_q = \\frac{\\sum_i p_i^q - 1}{(1-q)\\log b}
```
"""
tsallis_entropy(data, q::Real; base::Real=2, frequencies::Bool=true, species=nothing) =
    entropy(Tsallis(q; base), data; frequencies, species)

"""
    tsallis(data, q; base=2, frequencies=true)

Alias for [`tsallis_entropy`](@ref).
"""
tsallis(data, q::Real; base::Real=2, frequencies::Bool=true, species=nothing) =
    tsallis_entropy(data, q; base, frequencies, species)

"""
    tsallis_diversity(data, q; base=2, frequencies=true)

Return the effective diversity corresponding to Tsallis entropy of order `q`.

```math
{}^qD = \\left(1 + (1-q)(\\log b)T_q\\right)^{1/(1-q)}
```
"""
tsallis_diversity(data, q::Real; base::Real=2, frequencies::Bool=true, species=nothing) =
    effective_diversity(Tsallis(q; base), data; frequencies, species)

"""
    hill_number(data, q; frequencies=true)

Return the Hill number of order `q` for `data`.
"""
hill_number(data, q::Real; frequencies::Bool=true, species=nothing) =
    diversity(Hill(q), data; frequencies, species)

"""
    simpson_index(data; frequencies=true)

Return Simpson concentration, ``\\sum_i p_i^2``, for `data`.
"""
simpson_index(data; frequencies::Bool=true, species=nothing) =
    diversity(Simpson(), data; frequencies, species)

"""
    gini_simpson_index(data; frequencies=true)

Return Gini-Simpson diversity, ``1 - \\sum_i p_i^2``, for `data`.
"""
gini_simpson_index(data; frequencies::Bool=true, species=nothing) =
    diversity(GiniSimpson(), data; frequencies, species)

"""
    greenberg_diversity_index(data; frequencies=true)

Return Greenberg's linguistic diversity index, the probability that two
randomly selected individuals have different mother tongues.

This is equivalent to Gini-Simpson diversity, ``1 - \\sum_i p_i^2``.
"""
greenberg_diversity_index(data; frequencies::Bool=true, species=nothing) =
    diversity(GreenbergDiversityIndex(), data; frequencies, species)

"""
    linguistic_diversity_index(data; frequencies=true)

Return the linguistic diversity index (LDI), equivalent to
[`greenberg_diversity_index`](@ref).
"""
linguistic_diversity_index(data; frequencies::Bool=true, species=nothing) =
    diversity(LinguisticDiversityIndex(), data; frequencies, species)

"""
    index_of_linguistic_diversity(current, baseline; frequencies=true)

Return a Terralingua-style index of linguistic diversity (ILD) as the ratio of
the current linguistic diversity index to a baseline linguistic diversity index.

A value of `1` indicates no change relative to the baseline; values below `1`
indicate loss of linguistic diversity; values above `1` indicate an increase.
This helper compares two assemblages and is therefore not itself a
single-assemblage diversity index object.
"""
function index_of_linguistic_diversity(current, baseline; frequencies::Bool=true, species=nothing)
    baseline_value = linguistic_diversity_index(baseline; frequencies, species)
    baseline_value > 0 || throw(ArgumentError("baseline linguistic diversity must be positive"))
    return linguistic_diversity_index(current; frequencies, species) / baseline_value
end

"""
    inverse_simpson_index(data; frequencies=true)

Return inverse Simpson diversity, ``1 / \\sum_i p_i^2``, for `data`.
"""
inverse_simpson_index(data; frequencies::Bool=true, species=nothing) =
    diversity(InverseSimpson(), data; frequencies, species)
