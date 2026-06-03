# Getting Started

Load the package and create an assemblage as a dictionary from species to
abundance:

```julia
using DiversityIndices

assemblage = Dict(:oak => 12, :ash => 5, :elm => 3)
```

## Alpha Diversity

Convenience functions cover common indices:

```julia
richness(assemblage)
shannon_entropy(assemblage)
shannon_diversity(assemblage)
renyi_entropy(assemblage, 2)
tsallis_entropy(assemblage, 2)
simpson_index(assemblage)
gini_simpson_index(assemblage)
inverse_simpson_index(assemblage)
chao1(assemblage)
ace(assemblage)
sample_coverage(assemblage)
pielou_evenness(assemblage)
fisher_alpha(assemblage)
```

For notebook-style exploration, [`alpha_diversity`](@ref) collects the common
single-assemblage summaries into one named tuple:

```julia
alpha_diversity(assemblage)
```

Shannon, Renyi, and Tsallis entropy default to `base=2`, so entropy values are
reported in bits. Pass another `base` to the entropy index or convenience
function when you need a different logarithm base.

The generic API is useful when you want to choose an index explicitly. Use
`entropy` for entropy units and `diversity` for diversity / effective-diversity
units:

```julia
diversity(Richness(), assemblage)
entropy(Shannon(; base=2), assemblage)
entropy(Renyi(2), assemblage)
entropy(Tsallis(2), assemblage)
diversity(Shannon(; base=2), assemblage)
diversity(Renyi(2), assemblage)
diversity(Tsallis(2), assemblage)
effective_diversity(Shannon(; base=2), assemblage)
diversity(Hill(2), assemblage)
```

The main entropy and effective-diversity formulas are:

```math
H_b = -\sum_i p_i \log_b p_i
```

```math
H_q = \frac{1}{1-q}\log_b\left(\sum_i p_i^q\right)
```

```math
T_q = \frac{\sum_i p_i^q - 1}{(1-q)\log b}
```

```math
{}^qD = \left(\sum_i p_i^q\right)^{1/(1-q)}
```

## Input Conventions

Dictionaries are treated as `species => abundance` mappings.

```julia
proportions(Dict(:oak => 12, :ash => 5, :elm => 3))
```

Numeric vectors are treated as abundance vectors by default.

```julia
entropy(Shannon(), [12, 5, 3])
```

Non-numeric vectors are treated as raw observations.

```julia
observations = ["oak", "ash", "oak", "elm"]

counts(observations)
richness(observations)
```

Use `frequencies=false` when a numeric vector represents observations rather
than abundance values.

```julia
numeric_observations = [1, 2, 1, 3]

richness(numeric_observations; frequencies=false)
```

Community matrices are treated as samples in rows and taxa/categories in
columns. Alpha-diversity functions return one value per row.

```julia
community = [
    1 1 2 0 5
    3 0 1 1 0
]

richness(community)
entropy(Shannon(), community)
chao1(community)
ace(community)
sample_coverage(community)
alpha_diversity(community)
```

Tables.jl-compatible inputs, including DataFrames, use the same convention. Use
[`community_matrix`](@ref) to inspect the numeric matrix that will be analyzed.
By default numeric columns are used as species columns; pass `species` when the
table includes numeric metadata.

```julia
table = (
    site=["a", "b"],
    oak=[1, 3],
    ash=[1, 0],
    elm=[2, 1],
)

community_matrix(table; species=[:oak, :ash, :elm])
richness(table; species=[:oak, :ash, :elm])
shannon_entropy(table; species=[:oak, :ash, :elm])
```

## Pairwise Comparisons

Use [`similarity`](@ref) and [`dissimilarity`](@ref) with pairwise index types:

```julia
plot_a = Dict(:oak => 12, :ash => 5)
plot_b = Dict(:ash => 4, :elm => 7)

similarity(Jaccard(), plot_a, plot_b)
dissimilarity(Jaccard(), plot_a, plot_b)
similarity(SorensenDice(), plot_a, plot_b)
bray_curtis_dissimilarity(plot_a, plot_b)
```

When numeric vectors are used as abundance vectors for Bray-Curtis, positions
are treated as corresponding species and the vectors must have the same length.

```julia
left = [12, 5, 0]
right = [0, 4, 7]

dissimilarity(BrayCurtis(), left, right)
```

Passing a community matrix or Tables.jl-compatible table as the only data
argument returns a pairwise matrix across rows.

```julia
distance(BrayCurtis(), community)
jaccard_distance(community)
bray_curtis_distance(table; species=[:oak, :ash, :elm])
labeled_distance(BrayCurtis(), table; label=:site, species=[:oak, :ash, :elm])
```

Use the [similarity and dissimilarity catalog](similarity-index-catalog.md) for
index-selection notes, and the
[availability checklists](index-checklist.md) when comparing this package with
other biodiversity packages.

## Shannon Entropy Estimators

[`Shannon`](@ref) accepts an estimator object. The default is [`Plugin`](@ref),
the direct empirical estimator. Other estimators include [`MillerMadow`](@ref),
[`HausserStrimmer`](@ref), [`Basharin`](@ref), [`AddGamma`](@ref), and
[`ChaoShen`](@ref).

Use `support` when the finite support is known but some categories may have zero
observed counts. Pass an integer support size, or pass a collection of category
labels. Leave `support=nothing` when only the observed support is known, or when
using [`ChaoShen`](@ref) for possible unseen species.

```julia
entropy(Shannon(; estimator=Plugin()), assemblage)
entropy(Shannon(; estimator=MillerMadow()), assemblage)
entropy(Shannon(; estimator=AddGamma(1)), assemblage; support=10)
entropy(Shannon(; estimator=ChaoShen()), assemblage)

shannon_diversity(assemblage; estimator=MillerMadow())
shannon_diversity(assemblage; estimator=AddGamma(1), support=10)
```

## Generalized Entropies

[`Renyi`](@ref) and [`Tsallis`](@ref) provide entropy families indexed by order
``q``. The ``q = 1`` case is evaluated as Shannon entropy, and
[`effective_diversity`](@ref) returns the corresponding Hill number.

```julia
renyi_entropy(assemblage, 2)
renyi_diversity(assemblage, 2)

tsallis_entropy(assemblage, 2)
tsallis_diversity(assemblage, 2)
```

## Richness Estimation And Coverage

[`Chao1`](@ref), [`ACE`](@ref), and [`SampleCoverage`](@ref) provide
abundance-count estimators for possible unseen richness and sample completeness.

```julia
chao1(assemblage)
ace(assemblage)
ace(assemblage; threshold=10)
sample_coverage(assemblage)
```
