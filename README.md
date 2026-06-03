# DiversityAndDissimilarity.jl

[![CI](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/CI.yml)
[![Codecov](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/Codecov.yml/badge.svg)](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/Codecov.yml)
[![codecov](https://codecov.io/gh/mroughan/DiversityAndDissimilarity.jl/branch/main/graph/badge.svg?token=ZtXnAmM4b6)](https://codecov.io/gh/mroughan/DiversityAndDissimilarity.jl)
[![Aqua](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/Aqua.yml/badge.svg)](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/Aqua.yml)
[![JET](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/JET.yml/badge.svg)](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/JET.yml)
[![Documentation](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/mroughan/DiversityAndDissimilarity.jl/actions/workflows/documentation.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mroughan.github.io/DiversityAndDissimilarity.jl/dev)

`DiversityAndDissimilarity` provides ecology-style diversity indices for species
abundance data, linguistic-demography category counts, observation vectors,
and pairwise assemblage comparisons for dissimilarity or divergence.

The package is small and dependency-free. It is designed around dispatchable
index types such as `Richness()`, `Shannon()`, `Renyi()`, `Tsallis()`,
`Simpson()`, `Jaccard()`, `BrayCurtis()`, `Hellinger()`, and
`JensenShannon()`, plus convenience functions for common workflows.
Where other fields use the same mathematical object under a different name,
the package exposes those conventions explicitly; for example,
`LinguisticDiversityIndex()` and `GreenbergDiversityIndex()` are aliases by
interpretation for Gini-Simpson diversity.

If you are coming from R, the documentation includes a
[vegan migration guide](docs/src/vegan-migration.md) for the
[vegan](https://vegandevs.github.io/vegan/) package.
If you are choosing between packages, the documentation also includes
[availability checklists](docs/src/index-checklist.md) against common ecology
and biodiversity tools.

## Installation

From the Julia package prompt:

```julia
pkg> add DiversityAndDissimilarity
```

Or, while developing from this repository:

```julia
pkg> activate .
pkg> test
```

`DiversityAndDissimilarity` supports Julia 1.10 and newer 1.x releases.

## Quick Start

```julia
using DiversityAndDissimilarity

assemblage = Dict(:oak => 12, :ash => 5, :elm => 3)

richness(assemblage)
shannon_entropy(assemblage)
shannon_diversity(assemblage)
renyi_entropy(assemblage, 2)
tsallis_entropy(assemblage, 2)
simpson_index(assemblage)
chao1(assemblage)
sample_coverage(assemblage)
inverse_simpson_index(assemblage)
pielou_evenness(assemblage)
fisher_alpha(assemblage)
alpha_diversity(assemblage)
```

Shannon, Renyi, and Tsallis entropy default to `base=2`, so entropy values are
reported in bits. Pass another `base` to the entropy index or convenience
function when you need a different logarithm base.

For a more explicit style, pass an entropy index to `entropy` or an effective
diversity index to `diversity`:

```julia
diversity(Richness(), assemblage)
entropy(Shannon(; base=2), assemblage)
entropy(Renyi(2), assemblage)
entropy(Tsallis(2), assemblage)
diversity(Shannon(; base=2), assemblage)
diversity(Renyi(2), assemblage)
diversity(Tsallis(2), assemblage)
effective_diversity(Shannon(; base=2), assemblage)
effective_diversity(Renyi(2), assemblage)
diversity(Hill(1), assemblage)
```

## Data Inputs

Dictionaries are interpreted as `species => abundance`:

```julia
abundance = Dict("oak" => 12, "ash" => 5, "elm" => 3)
proportions(abundance)
```

Numeric vectors are interpreted as abundance vectors by default:

```julia
abundance = [12, 5, 3]
entropy(Shannon(), abundance)
```

Non-numeric vectors are interpreted as raw observations:

```julia
observations = ["oak", "ash", "oak", "elm"]
counts(observations)
richness(observations)
```

Pass `frequencies=false` to treat a numeric vector as raw observations rather
than abundance values:

```julia
numeric_observations = [1, 2, 1, 3]
richness(numeric_observations; frequencies=false)
```

Community matrices are interpreted as samples in rows and taxa/categories in
columns. Alpha-diversity functions return one value per row:

```julia
community = [
    1 1 2 0 5
    3 0 1 1 0
]

richness(community)
shannon_entropy(community)
chao1(community)
ace(community)
sample_coverage(community)
alpha_diversity(community)
```

Tables.jl-compatible inputs, including DataFrames, can be converted with
`community_matrix`. By default numeric columns are treated as species columns;
pass `species` explicitly when the table also contains numeric metadata.

```julia
community_matrix(table; species=[:oak, :ash, :elm])
richness(table; species=[:oak, :ash, :elm])
shannon_entropy(table; species=[:oak, :ash, :elm])
```

## Main API

The generic API is built around small index types:

- `entropy(index, data)` evaluates entropy-family indexes in entropy units.
- `diversity(index, data)` evaluates alpha diversity indexes; for entropy
  families it returns effective diversity.
- `effective_diversity(index, data)` returns the effective number of species
  where that transformation is standard.
- `similarity(index, left, right)` and `dissimilarity(index, left, right)`
  compare two assemblages.
- `index_metadata(index)`, `index_bounds(index)`, and trait helpers expose
  package conventions for generic workflows.

```julia
is_metric(Jaccard())
is_semimetric(BrayCurtis())
is_symmetric(KullbackLeibler())
index_bounds(JensenShannon())
```

The descriptor helpers include `is_finite`, `is_symmetric`, `is_nonnegative`,
`is_bounded`, `is_metric`, `is_triangular`, `is_pseudometric`,
`is_quasimetric`, `is_metametric`, `is_semimetric`, `is_premetric`,
`is_supermetric`, `is_similarity`, and `is_dissimilarity`. Unknown or
unencoded properties return `:unknown`.

Supported alpha diversity indices include:

- `Richness()`
- `Shannon(; base=2, estimator=Plugin())`
- `Renyi(q; base=2)`
- `Tsallis(q; base=2)`
- `Simpson()`
- `GiniSimpson()`
- `GreenbergDiversityIndex()`
- `LinguisticDiversityIndex()`
- `InverseSimpson()`
- `Hill(q)`
- `Chao1()`
- `ACE(; threshold=10)`
- `SampleCoverage()`
- `PielouEvenness()`
- `FisherAlpha()`

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

Supported pairwise indices include:

- `Jaccard()`
- `SorensenDice()`
- `Overlap()`
- `BrayCurtis()`
- `Ruzicka()`
- `TotalVariation()`
- `Manhattan()`
- `Euclidean()`
- `Canberra()`
- `Hellinger()`
- `Chord()`
- `Bhattacharyya()`
- `KullbackLeibler()`
- `ShannonDifference()`
- `JensenDifference()`
- `JensenShannon()`
- `MorisitaHorn()`

## Pairwise Comparisons

Incidence-style comparisons work with dictionaries or observation vectors:

```julia
plot_a = Dict(:oak => 12, :ash => 5)
plot_b = Dict(:ash => 4, :elm => 7)

similarity(Jaccard(), plot_a, plot_b)
jaccard_distance(plot_a, plot_b)
sorensen_dice_index(plot_a, plot_b)
overlap_similarity(plot_a, plot_b)
```

Abundance-style comparisons use aligned species abundances:

```julia
bray_curtis_dissimilarity(plot_a, plot_b)
ruzicka_similarity(plot_a, plot_b)
morisita_horn_similarity(plot_a, plot_b)
```

When passing numeric abundance vectors to `BrayCurtis()`, the vectors must have
the same length and positions are treated as shared species. Probability-style
distances normalize the aligned abundances first:

```julia
left = [12, 5, 0]
right = [0, 4, 7]

dissimilarity(BrayCurtis(), left, right)
hellinger_distance(left, right)
kullback_leibler_divergence(left, right)
shannon_difference(left, right)
jensen_difference(left, right)
jensen_shannon_similarity(left, right)
jensen_shannon_distance(left, right)
total_variation_distance(left, right)
```

As a rule of thumb, use Jaccard or Sorensen-Dice for presence/absence turnover,
Bray-Curtis for a practical ecological abundance default, Ruzicka when you want
a quantitative Jaccard, Hellinger or Chord for composition comparisons before
ordination, and Jensen-Shannon for a symmetric information-theoretic comparison.
Use Shannon difference when the question is whether two assemblages have
different entropy, and Jensen difference when you want the entropy gap between
the average composition and the average entropy of the two assemblages.
Kullback-Leibler divergence is available when an asymmetric
``D_{KL}(left \\Vert right)`` comparison is intended; it can be infinite when
the right assemblage has zero probability where the left assemblage is
positive.

Low-sample corrections for KL and Jensen-Shannon-style divergences use the same
estimator objects as Shannon entropy:

```julia
kullback_leibler_divergence(left, right; estimator=MillerMadow())
kullback_leibler_divergence(left, right; estimator=AddGamma(1))    # Laplace
kullback_leibler_divergence(left, right; estimator=AddGamma(0.5))  # Jeffreys
kullback_leibler_divergence(left, right; estimator=HausserStrimmer())
kullback_leibler_divergence(left, right; estimator=ChaoShen())

jensen_shannon_divergence(left, right; estimator=MillerMadow())
jensen_shannon_distance(left, right; estimator=AddGamma(0.5), support=10)
```

`AddGamma` and `HausserStrimmer` use the aligned observed support unless a
larger known `support` is supplied. `ChaoShen()` applies a Good-Turing-style
unseen-mass correction.

Passing a community matrix or Tables.jl-compatible table as the only data
argument returns a pairwise matrix across rows:

```julia
distance(BrayCurtis(), community)
distance(Hellinger(), community)
jensen_shannon_distance(community)
jaccard_distance(community)
bray_curtis_distance(table; species=[:oak, :ash, :elm])
labeled_distance(BrayCurtis(), table; label=:site, species=[:oak, :ash, :elm])
```

For scaling, alpha-diversity summaries are generally linear in observed support
or in the number of entries in a community matrix. Full pairwise distance
matrices scale quadratically in the number of samples/sites because they return
an `M x M` result. Known-support estimators may also scale with the supplied
support size, and bootstrap intervals scale roughly linearly with `nboot`.

## Shannon Entropy Estimators

Shannon entropy estimation is represented by dispatchable estimator types. The
default is `Plugin()`, the direct empirical estimator. Other estimators include
`MillerMadow()`, `HausserStrimmer()`, `Basharin()`, `AddGamma(gamma)`, and
`ChaoShen()`.

Use `support` when the finite support is known but some categories may have zero
observed counts. Pass an integer support size, or pass a collection of category
labels. Leave `support=nothing` when only the observed support is known, or when
using `ChaoShen()` for possible unseen species.

```julia
entropy(Shannon(; estimator=Plugin()), assemblage)
entropy(Shannon(; estimator=MillerMadow()), assemblage)
entropy(Shannon(; estimator=AddGamma(1)), assemblage; support=10)
entropy(Shannon(; estimator=ChaoShen()), assemblage)

shannon_diversity(assemblage; estimator=MillerMadow())
shannon_diversity(assemblage; estimator=AddGamma(1), support=10)
```

Analytic or approximate Shannon uncertainty is available for selected
estimators, and bootstrap/jackknife helpers are available for integer count
data:

```julia
shannon_variance(assemblage; estimator=Basharin(), support=10)
shannon_confint(assemblage; estimator=ChaoShen())
bootstrap(Shannon(; estimator=AddGamma(1)), assemblage; support=10)
jackknife(Shannon(; estimator=MillerMadow()), assemblage)
```

## Generalized Entropies

Renyi and Tsallis entropies are available as dispatchable index types and
convenience functions. The ``q = 1`` case is evaluated as Shannon entropy, and
`effective_diversity` returns the corresponding Hill number.

```julia
renyi_entropy(assemblage, 2)
renyi_diversity(assemblage, 2)

tsallis_entropy(assemblage, 2)
tsallis_diversity(assemblage, 2)
```

## Richness Estimation And Coverage

Chao1, ACE, and Good-Turing sample coverage are available for abundance count
data:

```julia
chao1(assemblage)
ace(assemblage)
ace(assemblage; threshold=10)
sample_coverage(assemblage)
```

These functions also work row-wise on community matrices.

## Documentation

[Documenter.jl](https://documenter.juliadocs.org/stable/) documentation lives
in `docs/`.

```julia
julia --project=docs docs/make.jl
```

The generated site is written to `docs/build/`.

Useful entry points include the [getting-started guide](docs/src/getting-started.md),
the [vegan migration guide](docs/src/vegan-migration.md), the
[index availability checklists](docs/src/index-checklist.md), and the
[benchmark notes](docs/src/benchmarks.md).

## Benchmarks

Reproducible benchmark scripts live in `benchmark/`. They generate a
deterministic simulated community matrix at runtime, so no large dataset is
committed and Git LFS is not needed.

```bash
julia --project=. benchmark/julia_benchmark.jl
python3 benchmark/python_benchmark.py
Rscript benchmark/r_vegan_benchmark.R
benchmark/run_report.sh
```

The Python benchmark requires [NumPy](https://numpy.org/) and uses
[SciPy](https://scipy.org/) for Bray-Curtis when available. The R benchmark
requires [vegan](https://vegandevs.github.io/vegan/).

## Quality Checks

The ordinary package tests are run with:

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

Additional package-quality checks use
[Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) and
[JET.jl](https://github.com/aviatesk/JET.jl) from a dedicated test environment:

```julia
julia --project=test/quality -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=test/quality test/quality/runtests.jl
```

GitHub Actions runs the standard Julia compatibility matrix separately from
the Aqua and JET quality workflows, which run only on the latest stable Julia
release.

## Disclosure

This package was developed with assistance from OpenAI Codex, an AI coding
assistant based on GPT-5. Code design decisions were human mediated, and the
resulting code was manually reviewed.
