# Diversity Indices

Diversity indices summarize one assemblage at a time. In
`DiversityAndDissimilarity.jl`, use [`entropy`](@ref) when you want entropy
units and [`diversity`](@ref) or [`effective_diversity`](@ref) when you want
diversity-scale values.

```jldoctest diversitypage
julia> using DiversityAndDissimilarity

julia> assemblage = Dict(:oak => 12, :ash => 5, :elm => 3);

julia> diversity(Richness(), assemblage)
3

julia> entropy(Shannon(), assemblage)
1.3527241956246545

julia> diversity(Shannon(), assemblage)
2.5539392274300625
```

## Available Methods

| Family | Types and helpers | Main options | Typical use |
|---|---|---|---|
| Richness | [`Richness`](@ref), [`richness`](@ref) | `frequencies` | Observed positive-abundance taxa/categories. |
| Shannon entropy | [`Shannon`](@ref), [`shannon_entropy`](@ref), [`shannon_diversity`](@ref) | `base`, `estimator`, `support` | Entropy and effective diversity; low-count corrections. |
| Generalized entropy | [`Renyi`](@ref), [`Tsallis`](@ref), [`renyi_entropy`](@ref), [`tsallis_entropy`](@ref) | order `q`, `base` | Diversity profiles and Hill-number relationships. |
| Simpson family | [`Simpson`](@ref), [`GiniSimpson`](@ref), [`InverseSimpson`](@ref) | `frequencies` | Dominance, Gini-Simpson diversity, and inverse Simpson diversity. |
| Hill numbers | [`Hill`](@ref), [`hill_number`](@ref), [`effective_diversity`](@ref) | order `q` | Effective number of species/categories. |
| Richness estimators | [`Chao1`](@ref), [`ACE`](@ref), [`chao1`](@ref), [`ace`](@ref) | ACE `threshold` | Unseen-richness estimation from abundance frequency counts. |
| Coverage | [`SampleCoverage`](@ref), [`sample_coverage`](@ref) | `frequencies` | Good-Turing style sample completeness. |
| Evenness | [`PielouEvenness`](@ref), [`pielou_evenness`](@ref) | `base` | Shannon evenness relative to observed richness. |
| Fisher log-series | [`FisherAlpha`](@ref), [`fisher_alpha`](@ref) | `frequencies` | Fisher's alpha diversity parameter. |
| Linguistic diversity | [`GreenbergDiversityIndex`](@ref), [`LinguisticDiversityIndex`](@ref), [`linguistic_diversity_index`](@ref) | `frequencies` | Gini-Simpson interpreted as probability of different mother tongues. |

## Inputs

All methods accept the same data shapes: dictionaries, numeric vectors,
observation vectors, community matrices, and Tables.jl-compatible tables.
See [Data Input Formats](data-input.md) for the complete reference.

```jldoctest diversitypage
julia> richness(["oak", "ash", "oak", "elm"])   # observation vector
3

julia> richness([1, 2, 1, 3]; frequencies=false)   # numeric obs. vector
3

julia> community = [1 1 2 0 5; 3 0 1 1 0];

julia> all(shannon_entropy(community) .≈ [1.6577427265048888, 1.3709505944546687])
true
```

## Entropy And Effective Diversity

[`Shannon`](@ref), [`Renyi`](@ref), and [`Tsallis`](@ref) default to `base=2`,
so entropy is reported in bits unless you choose another base.

```math
H_b = -\sum_i p_i \log_b p_i
```

```math
H_q = \frac{1}{1-q}\log_b\left(\sum_i p_i^q\right)
```

```math
T_q = \frac{\sum_i p_i^q - 1}{(1-q)\log b}
```

The order `q = 1` case for Renyi and Tsallis is evaluated as Shannon entropy.
Effective diversity converts entropy-family quantities to Hill-number scale:

```math
{}^qD = \left(\sum_i p_i^q\right)^{1/(1-q)}, \qquad {}^1D = b^H.
```

```jldoctest diversitypage
julia> entropy(Renyi(2), assemblage)
1.168122758808327

julia> effective_diversity(Renyi(2), assemblage)
2.247191011235955

julia> diversity(Hill(2), assemblage)
2.247191011235955
```

## Shannon Estimators And Low-Count Options

[`Shannon`](@ref) accepts estimator objects:

```jldoctest diversitypage
julia> entropy(Shannon(; estimator=Plugin()), assemblage)
1.3527241956246545

julia> entropy(Shannon(; estimator=MillerMadow()), assemblage)
1.4248589476691027

julia> entropy(Shannon(; estimator=AddGamma(1)), assemblage; support=5)
1.7792365361682794
```

Use [`Plugin`](@ref) for the empirical estimate. Use [`MillerMadow`](@ref) for
a simple observed-support bias correction. Use [`HausserStrimmer`](@ref),
[`Basharin`](@ref), or [`AddGamma`](@ref) when a finite support is known. Use
[`ChaoShen`](@ref) when support is unknown and unseen categories are plausible.

The `support` keyword can be an integer support size or a collection of known
category labels:

```jldoctest diversitypage
julia> known_species = [:oak, :ash, :elm, :pine, :birch];

julia> entropy(Shannon(; estimator=AddGamma(0.5)), assemblage; support=known_species)
1.6295944962456386
```

For uncertainty, use analytic helpers where available or resampling:

```julia
entropy_variance(Shannon(; estimator=Basharin()), assemblage; support=5)
entropy_confint(Shannon(; estimator=ChaoShen()), assemblage)
bootstrap(Shannon(; estimator=AddGamma(1)), assemblage; support=5)
jackknife(Shannon(; estimator=MillerMadow()), assemblage)
```

## Simpson, Linguistic Diversity, Richness, And Coverage

The Simpson family uses concentration

```math
D = \sum_i p_i^2.
```

[`Simpson`](@ref) returns `D`, [`GiniSimpson`](@ref) returns `1 - D`, and
[`InverseSimpson`](@ref) returns `1 / D`.

[`GreenbergDiversityIndex`](@ref) and [`LinguisticDiversityIndex`](@ref) are
the same formula as Gini-Simpson but interpreted as the probability that two
randomly selected people have different mother tongues.

Richness and coverage estimators use abundance-frequency counts:

```jldoctest diversitypage
julia> sparse = [5, 1, 1, 0, 0];

julia> richness(sparse)
3

julia> chao1(sparse)
4.0

julia> sample_coverage(sparse)
0.7142857142857143
```

## Availability Checklist

Legend: `[x]` is available directly; `[~]` is available indirectly or with a
different convention; `[ ]` is not documented as available in the checked
source. Last checked: 2026-05-13.

| Index or feature | DiversityAndDissimilarity.jl | Diversity.jl | vegan | iNEXT | scikit-bio | EcoPy | Microbiome.jl | SciPy | SpadeR | entropart | Notes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Observed richness / number of taxa | [x] | [x] | [x] | [x] | [x] | [x] | [ ] | [ ] | [x] | [x] | `Richness()` / `richness`; vegan `specnumber`; scikit-bio `sobs` / `observed_features`; EcoPy `spRich`. |
| Chao / ACE richness estimators | [x] | [ ] | [x] | [x] | [x] | [ ] | [ ] | [ ] | [x] | [~] | `Chao1()` / `chao1` and `ACE()` / `ace`; vegan `estimateR` / `specpool`; iNEXT `ChaoRichness`; scikit-bio `chao1` and `ace`. |
| Rarefied richness / rarefaction curves | [ ] | [ ] | [x] | [x] | [~] | [x] | [ ] | [ ] | [~] | [ ] | Outside this package's core scope; vegan and iNEXT are stronger rarefaction tools. |
| Shannon entropy, `H` | [x] | [~] | [x] | [~] | [x] | [x] | [x] | [ ] | [x] | [x] | Defaults to `base=2` here; vegan defaults to natural logs. |
| Shannon effective diversity | [x] | [x] | [~] | [x] | [x] | [x] | [ ] | [ ] | [x] | [x] | `effective_diversity(Shannon())`; vegan derives it with `exp(diversity(..., "shannon"))`. |
| Shannon entropy estimation | [x] | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [x] | [x] | `Plugin`, `MillerMadow`, `HausserStrimmer`, `Basharin`, `AddGamma`, and `ChaoShen` are available here. |
| Sample coverage | [x] | [ ] | [~] | [x] | [x] | [ ] | [ ] | [ ] | [x] | [x] | `SampleCoverage()` / `sample_coverage`; related coverage estimates appear in iNEXT, scikit-bio, SpadeR, and entropart. |
| Simpson concentration, `sum(p_i^2)` | [x] | [~] | [~] | [~] | [~] | [x] | [ ] | [ ] | [x] | [x] | This package's `Simpson()` returns concentration; vegan `index="simpson"` returns `1-D`. |
| Gini-Simpson, `1 - sum(p_i^2)` | [x] | [~] | [x] | [~] | [x] | [x] | [x] | [ ] | [x] | [~] | `GiniSimpson()`; vegan `diversity(..., "simpson")`. |
| Greenberg / linguistic diversity index | [x] | [~] | [~] | [~] | [~] | [~] | [~] | [ ] | [~] | [~] | Same formula as Gini-Simpson with linguistic-demography interpretation. |
| Inverse Simpson | [x] | [x] | [x] | [x] | [x] | [x] | [ ] | [ ] | [x] | [x] | `InverseSimpson()`; vegan `index="invsimpson"`; scikit-bio `inv_simpson`. |
| Hill number, general order `q` | [x] | [x] | [~] | [x] | [x] | [~] | [ ] | [ ] | [x] | [x] | `Hill(q)` and `hill_number`; Diversity.jl is stronger for partitioned Hill/Jost diversity. |
| Renyi entropy / diversity profile | [x] | [~] | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [~] | [x] | `Renyi(q; base=2)` / `renyi_entropy`; vegan `renyi`; scikit-bio `renyi`. |
| Tsallis entropy | [x] | [~] | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [x] | `Tsallis(q; base=2)` / `tsallis_entropy`; vegan `tsallis`; scikit-bio `tsallis`. |
| Fisher alpha | [x] | [ ] | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | `FisherAlpha()` / `fisher_alpha`; vegan `fisher.alpha`; scikit-bio `fisher_alpha`. |
| Evenness / equitability | [x] | [ ] | [~] | [ ] | [x] | [x] | [ ] | [ ] | [ ] | [ ] | `PielouEvenness()` / `pielou_evenness`; scikit-bio and EcoPy provide additional evenness metrics. |
| Dominance / Berger-Parker style metrics | [ ] | [ ] | [~] | [ ] | [x] | [x] | [ ] | [ ] | [ ] | [ ] | Simpson concentration is available here, but broader dominance families are not. |
| UniFrac / Faith phylogenetic diversity | [ ] | [x] | [~] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [x] | Out of scope here; Diversity.jl and scikit-bio are stronger phylogenetic options. |
| Rao quadratic entropy | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [~] | Specialist phylogenetic/functional diversity packages are better suited. |
| Alpha/beta/gamma diversity partitioning | [ ] | [x] | [x] | [~] | [~] | [x] | [ ] | [ ] | [~] | [x] | Not a core goal of this package. |
| Functional diversity | [ ] | [~] | [~] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [x] | Specialist packages such as BAT, hillR, hilldiv, and entropart are stronger. |

## Reference

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
