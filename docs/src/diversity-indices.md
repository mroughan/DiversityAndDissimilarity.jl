# Diversity Indices

Diversity and entropy indices summarize a single assemblage. In
`DiversityAndDissimilarity`, use [`entropy`](@ref) for entropy units and
[`diversity`](@ref) for diversity / effective-diversity units:

```julia
diversity(Richness(), assemblage)
entropy(Shannon(), assemblage)
entropy(Renyi(2), assemblage)
entropy(Tsallis(2), assemblage)
diversity(Shannon(), assemblage)
diversity(Renyi(2), assemblage)
diversity(Tsallis(2), assemblage)
diversity(Simpson(), assemblage)
```

The input is converted to positive relative abundances, ``p_i``. Dictionaries
are interpreted as `species => abundance`, numeric vectors are interpreted as
abundance vectors by default, and non-numeric vectors are interpreted as raw
observations. Community matrices are interpreted as samples in rows and
taxa/categories in columns; alpha-diversity operations return one value per row.
Tables.jl-compatible inputs, including DataFrames, can be used directly; pass
`species` to choose the species columns when needed.

```julia
community_matrix(table; species=[:oak, :ash, :elm])
richness(table; species=[:oak, :ash, :elm])
entropy(Shannon(), table; species=[:oak, :ash, :elm])
```

For one assemblage, most alpha-diversity calculations are linear in observed
support, ``O(S)``. For a dense community matrix with ``M`` samples and ``P``
species columns, row-wise alpha diversity is ``O(MP)``. See
[Scaling And Performance](@ref) for more detail on support size, community
matrix size, and pairwise distance matrices.

## Entropy Families

[`Shannon`](@ref), [`Renyi`](@ref), and [`Tsallis`](@ref) are entropy-based
indices. Use [`entropy`](@ref) to return entropy values. They default to
`base=2`, so entropy is measured in bits unless you choose another base.

```math
H_b = -\\sum_i p_i \\log_b p_i
```

```math
H_q = \\frac{1}{1-q}\\log_b\\left(\\sum_i p_i^q\\right)
```

```math
T_q = \\frac{\\sum_i p_i^q - 1}{(1-q)\\log b}
```

The order ``q = 1`` case for Renyi and Tsallis entropy is evaluated as Shannon
entropy.

!!! note "Tsallis base scaling"
    The ``\\log b`` denominator in the Tsallis formula above scales the result so
    that ``T_1`` equals Shannon entropy in the same base. This differs from the
    standard textbook definition ``(1 - \\sum_i p_i^q)/(q-1)``, which has no
    logarithm base factor. Values will differ from packages using the standard
    definition by a factor of ``\\log b``.

## Shannon Entropy Estimators

[`Shannon`](@ref) estimates entropy using an estimator object:

```julia
entropy(Shannon(; estimator=Plugin()), assemblage)
entropy(Shannon(; estimator=MillerMadow()), assemblage)
entropy(Shannon(; estimator=HausserStrimmer()), assemblage; support=10)
entropy(Shannon(; estimator=Basharin()), assemblage; support=10)
entropy(Shannon(; estimator=AddGamma(1)), assemblage; support=10)
entropy(Shannon(; estimator=ChaoShen()), assemblage)

shannon_diversity(assemblage; estimator=MillerMadow())
shannon_diversity(assemblage; estimator=AddGamma(1), support=10)
```

When the finite support is known, pass it as `support`. This may be an integer
support size or a collection of category labels. When the support is unknown
and unseen categories may exist, leave `support=nothing`; [`ChaoShen`](@ref) is
designed for that setting.

Known-support estimators can scale with the supplied support size ``K``, not
only with observed support ``S``. This matters when the possible species list is
very large relative to the species observed in a sample.

See [Choosing An Estimator](@ref) for a flowchart-based guide to selecting and
testing estimators on the same assemblage.

## Uncertainty

For [`Plugin`](@ref), [`MillerMadow`](@ref), [`Basharin`](@ref), and
[`ChaoShen`](@ref), use [`entropy_variance`](@ref) or
[`entropy_confint`](@ref) for analytic or approximate Shannon entropy
uncertainty estimates:

```julia
entropy_variance(Shannon(; estimator=Basharin()), assemblage; support=10)
entropy_confint(Shannon(; estimator=ChaoShen()), assemblage)

shannon_variance(assemblage; estimator=Basharin(), support=10)
shannon_confint(assemblage; estimator=ChaoShen())
```

For estimators without an analytic variance in this package, use resampling:

```julia
bootstrap(Shannon(; estimator=AddGamma(1)), assemblage; support=10)
jackknife(Shannon(; estimator=MillerMadow()), assemblage)
bootstrap(Shannon(), assemblage; quantity=:diversity)
```

Bootstrap and jackknife helpers require integer count data. They work row-wise
for community matrices and Tables.jl-compatible inputs.

Resampling repeats the underlying estimator many times. Bootstrap cost grows
roughly linearly with `nboot`, while jackknife cost grows with the sample size
or total integer count represented by the assemblage.

## Effective Diversity

[`effective_diversity`](@ref) converts supported indices into effective species
numbers. For entropy-family indices this returns the corresponding Hill number:

```math
{}^qD = \\left(\\sum_i p_i^q\\right)^{1/(1-q)}.
```

For Shannon entropy this is:

```math
{}^1D = b^H.
```

## Simpson Family

The Simpson-family indices are different transformations of the same
concentration quantity:

```math
D = \\sum_i p_i^2.
```

- [`Simpson`](@ref) returns ``D``.
- [`GiniSimpson`](@ref) returns ``1-D``.
- [`InverseSimpson`](@ref) returns ``1/D``.

## Linguistic Diversity

[`GreenbergDiversityIndex`](@ref) and [`LinguisticDiversityIndex`](@ref) are
linguistic-demography names for the same quantity as [`GiniSimpson`](@ref):

```math
LDI = 1 - \sum_i p_i^2.
```

Interpreted linguistically, ``p_i`` is the proportion of people with mother
tongue ``i``. The value is the probability that two randomly selected people
have different mother tongues. This matches Greenberg's linguistic diversity
index and is included because it has the same mathematical structure as an
ecological diversity index.

```julia
greenberg_diversity_index(language_counts)
linguistic_diversity_index(language_counts)
diversity(LinguisticDiversityIndex(), language_counts)
```

The index does not account for second-language use, language vitality, or how
different the languages are from each other. It only uses the relative
frequencies of language categories.

The temporal Index of Linguistic Diversity can be computed as a ratio relative
to a baseline:

```julia
index_of_linguistic_diversity(current_counts, baseline_counts)
```

A value below `1` indicates loss relative to the baseline; values above `1`
indicate an increase. This helper compares two assemblages, so it is not
itself a single-assemblage `DiversityIndex`.

## Richness Estimation And Coverage

[`Chao1`](@ref), [`ACE`](@ref), and [`SampleCoverage`](@ref) use abundance-count
frequencies to estimate unseen richness or sample completeness.

For observed richness ``S_{obs}``, singleton count ``f_1``, doubleton count
``f_2``, and sample size ``n``:

```math
\\hat S_{Chao1} = S_{obs} + \\frac{f_1(f_1 - 1)}{2(f_2 + 1)}
```

```math
\\hat C = 1 - \\frac{f_1}{n}
```

ACE splits taxa into rare and abundant groups using a count threshold, defaulting
to `threshold=10`.

```julia
chao1(assemblage)
ace(assemblage)
ace(assemblage; threshold=10)
sample_coverage(assemblage)
```

## Evenness And Fisher Alpha

[`PielouEvenness`](@ref) reports Shannon evenness relative to the maximum
entropy for the observed richness:

```math
J = \\frac{H}{\\log S}.
```

[`FisherAlpha`](@ref) solves Fisher's log-series relationship between observed
richness ``S`` and total abundance ``n``:

```math
S = \\alpha \\log\\left(1 + \\frac{n}{\\alpha}\\right).
```

```julia
pielou_evenness(assemblage)
fisher_alpha(assemblage)
```

## Convenience Functions

Convenience functions call the same generic methods:

```julia
richness(assemblage)
shannon_entropy(assemblage)
shannon(assemblage)
renyi_entropy(assemblage, 2)
renyi(assemblage, 2)
tsallis_entropy(assemblage, 2)
tsallis(assemblage, 2)
hill_number(assemblage, 2)
simpson_index(assemblage)
inverse_simpson_index(assemblage)
pielou_evenness(assemblage)
fisher_alpha(assemblage)
```

```@docs
ShannonEstimator
Plugin
MillerMadow
HausserStrimmer
Basharin
AddGamma
ChaoShen
Richness
Shannon
Renyi
Tsallis
Simpson
GiniSimpson
GreenbergDiversityIndex
LinguisticDiversityIndex
InverseSimpson
Hill
Chao1
ACE
SampleCoverage
PielouEvenness
FisherAlpha
entropy
entropy_variance
entropy_confint
diversity
effective_diversity
richness
shannon
shannon_entropy
shannon_variance
shannon_confint
shannon_diversity
bootstrap
jackknife
renyi
renyi_entropy
renyi_diversity
tsallis
tsallis_entropy
tsallis_diversity
hill_number
chao1
ace
sample_coverage
pielou_evenness
fisher_alpha
simpson_index
gini_simpson_index
greenberg_diversity_index
linguistic_diversity_index
index_of_linguistic_diversity
inverse_simpson_index
```
