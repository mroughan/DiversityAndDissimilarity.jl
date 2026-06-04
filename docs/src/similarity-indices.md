# Dissimilarity Indices

Pairwise indices compare two assemblages. Some are naturally similarities,
some are dissimilarities or distances, and some are information divergences.
`DiversityAndDissimilarity.jl` exposes the convention explicitly through
[`similarity`](@ref), [`dissimilarity`](@ref), and [`distance`](@ref).

```jldoctest dissimilaritypage
julia> using DiversityAndDissimilarity

julia> left = Dict(:oak => 12, :ash => 5);

julia> right = Dict(:ash => 4, :elm => 7);

julia> similarity(Jaccard(), left, right)
0.3333333333333333

julia> dissimilarity(BrayCurtis(), left, right)
0.7142857142857143
```

## Available Methods

| Family | Types and helpers | Main options | Typical use |
|---|---|---|---|
| Incidence | [`Jaccard`](@ref), [`SorensenDice`](@ref), [`Overlap`](@ref) | `frequencies` | Presence/absence overlap and turnover. |
| Abundance | [`BrayCurtis`](@ref), [`Ruzicka`](@ref), [`Canberra`](@ref), [`MorisitaHorn`](@ref) | `frequencies` | Ecological abundance comparisons with aligned taxa. |
| Probability distances | [`TotalVariation`](@ref), [`Manhattan`](@ref), [`Euclidean`](@ref), [`Hellinger`](@ref), [`Chord`](@ref), [`Bhattacharyya`](@ref) | `frequencies` | Comparisons of normalized composition. |
| Information divergences | [`KullbackLeibler`](@ref), [`ShannonDifference`](@ref), [`JensenDifference`](@ref), [`JensenShannon`](@ref) | `base`, `estimator`, `support`, `distance` | Directional KL, entropy differences, and Jensen-Shannon divergence/distance. |
| Matrix helpers | [`labeled_distance`](@ref), [`labeled_dissimilarity`](@ref), [`labeled_similarity`](@ref) | `labels`, `label`, `species` | Pairwise matrices with sample labels. |

## Choosing A Measure

- Use [`Jaccard`](@ref) for simple presence/absence turnover.
- Use [`SorensenDice`](@ref) when shared presences should receive more weight
  than in Jaccard.
- Use [`Overlap`](@ref) when asking whether a smaller assemblage is nested in a
  larger one.
- Use [`BrayCurtis`](@ref) as a robust default for ecological abundance data.
- Use [`Ruzicka`](@ref) when you want an abundance analogue of Jaccard.
- Use [`Hellinger`](@ref) or [`Chord`](@ref) before Euclidean-style workflows
  on transformed relative abundances.
- Use [`TotalVariation`](@ref) for a bounded probability-composition
  difference with direct mass interpretation.
- Use [`JensenShannon`](@ref) for a symmetric information-theoretic metric.
- Use [`KullbackLeibler`](@ref) only when the direction
  `left || right` is scientifically meaningful.

## Incidence Examples

Incidence comparisons use species presence or absence:

```math
J(A,B) = \frac{|A \cap B|}{|A \cup B|}, \qquad
S(A,B) = \frac{2|A \cap B|}{|A| + |B|}.
```

```jldoctest dissimilaritypage
julia> a = [1, 1, 0, 1];

julia> b = [1, 0, 1, 1];

julia> similarity(Jaccard(), a, b)
0.5

julia> dissimilarity(Jaccard(), a, b)
0.5

julia> similarity(SorensenDice(), a, b)
0.6666666666666666
```

## Abundance Examples

Abundance comparisons use aligned abundance vectors. For dictionaries, taxa are
aligned by key. For numeric vectors, positions are corresponding taxa.

```math
BC(x,y) = \frac{\sum_i |x_i-y_i|}{\sum_i (x_i+y_i)}.
```

```jldoctest dissimilaritypage
julia> x = [1, 2, 3];

julia> y = [2, 2, 0];

julia> bray_curtis_dissimilarity(x, y)
0.4

julia> ruzicka_similarity(x, y)
0.42857142857142855

julia> canberra_distance(x, y)
0.4444444444444444
```

[`MorisitaHorn`](@ref) is an abundance-overlap similarity dominated by shared
common taxa:

```julia
morisita_horn_similarity(x, y)
morisita_horn_distance(x, y)
```

## Probability And Information Examples

Probability comparisons normalize aligned abundances to `p` and `q`.

```math
TV(p,q) = \frac{1}{2}\sum_i |p_i-q_i|.
```

```math
D_{KL}(p \Vert q) = \sum_i p_i \log_b\frac{p_i}{q_i}.
```

```math
JS(p,q) = \frac{1}{2}D_{KL}(p \Vert m) +
          \frac{1}{2}D_{KL}(q \Vert m), \qquad m = \frac{p+q}{2}.
```

```jldoctest dissimilaritypage
julia> total_variation_distance([1, 0, 0], [0, 1, 0])
1.0

julia> hellinger_distance([1, 0, 0], [0, 1, 0])
1.0

julia> jensen_shannon_distance([1, 0], [0, 1])
1.0

julia> dissimilarity(KullbackLeibler(), [1, 0], [0, 1])
Inf
```

[`KullbackLeibler`](@ref) is asymmetric, so pairwise matrices are not forced to
be symmetric. [`JensenShannon`](@ref) returns the square-root distance by
default; use `JensenShannon(; distance=false)` or [`JensenDifference`](@ref)
for the raw divergence.

## Low-Sample Divergence Corrections

KL, Jensen difference, and Jensen-Shannon divergence/distance accept the same
estimator objects used for Shannon entropy:

```julia
kullback_leibler_divergence(left, right; estimator=MillerMadow())
kullback_leibler_divergence(left, right; estimator=AddGamma(1))    # Laplace
kullback_leibler_divergence(left, right; estimator=AddGamma(0.5))  # Jeffreys
kullback_leibler_divergence(left, right; estimator=HausserStrimmer())
kullback_leibler_divergence(left, right; estimator=ChaoShen())

jensen_shannon_divergence(left, right; estimator=AddGamma(0.5), support=10)
jensen_shannon_distance(left, right; estimator=ChaoShen())
```

Use `support` when the finite category universe is known. `AddGamma(1)` is
Laplace smoothing; `AddGamma(0.5)` is Jeffreys smoothing. `ChaoShen()` applies
sample-coverage logic for unseen mass.

## Community Matrices

Passing a community matrix computes all pairwise comparisons across rows:

```jldoctest dissimilaritypage
julia> community = [
           1 1 2 0 5
           3 0 1 1 0
       ];

julia> distance(BrayCurtis(), community)
2×2 Matrix{Float64}:
 0.0       0.714286
 0.714286  0.0

julia> labeled_distance(BrayCurtis(), community; labels=["plot-a", "plot-b"]).labels
2-element Vector{String}:
 "plot-a"
 "plot-b"
```

For Tables.jl-compatible inputs, use `species` to select taxa columns and
`label` to carry sample identifiers.

## Availability Checklist

Legend: `[x]` is available directly; `[~]` is available indirectly or with a
different convention; `[ ]` is not documented as available in the checked
source. Last checked: 2026-05-15.

| Index or feature | DiversityAndDissimilarity.jl | Diversity.jl | vegan | iNEXT | scikit-bio | EcoPy | Microbiome.jl | SciPy | SpadeR | entropart / hill packages | Notes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Jaccard incidence similarity/distance | [x] | [ ] | [x] | [ ] | [~] | [x] | [x] | [x] | [x] | [ ] | `Jaccard()` / `jaccard_index`; vegan `vegdist(method="jaccard")`; SciPy boolean `jaccard`. |
| Sorensen-Dice incidence similarity/distance | [x] | [ ] | [~] | [ ] | [~] | [x] | [ ] | [x] | [x] | [ ] | `SorensenDice()`; SciPy boolean `dice`; related vegan forms exist through binary transformations. |
| Overlap / Szymkiewicz-Simpson | [x] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [ ] | `Overlap()` / `overlap_similarity`; nestedness-sensitive incidence overlap. |
| Simple matching / Sokal-Michener | [ ] | [ ] | [ ] | [ ] | [~] | [x] | [ ] | [x] | [ ] | [ ] | Requires meaningful shared absences and a fixed species universe. |
| Russell-Rao | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [x] | [ ] | [ ] | Shared-absence-sensitive boolean coefficient. |
| Kulczynski / Mountford / Raup-Crick incidence | [ ] | [ ] | [x] | [ ] | [~] | [ ] | [ ] | [ ] | [~] | [ ] | Available in vegan or specialist packages; conventions vary. |
| Chao-Jaccard / Chao-Sorensen | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [x] | [~] | Adjusts overlap for unseen shared species. |
| Bray-Curtis dissimilarity | [x] | [ ] | [x] | [ ] | [~] | [x] | [x] | [x] | [~] | [ ] | `BrayCurtis()`; vegan `vegdist(method="bray")`; SciPy `braycurtis`. |
| Quantitative Sorensen / Bray-Curtis similarity | [~] | [ ] | [~] | [ ] | [~] | [x] | [~] | [~] | [~] | [ ] | Derivable as `1 - bray_curtis_distance`; naming conventions vary. |
| Ruzicka / quantitative Jaccard | [x] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [ ] | `Ruzicka()` / `ruzicka_similarity`; abundance version of Jaccard. |
| Percentage similarity / Renkonen | [ ] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | Probability overlap `sum(min(p,q))`; complement is total variation. |
| Total variation distance | [x] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [~] | [ ] | [ ] | `TotalVariation()` / `total_variation_distance`; equals Bray-Curtis for normalized probabilities. |
| Manhattan / L1 distance | [x] | [ ] | [x] | [ ] | [~] | [x] | [ ] | [x] | [ ] | [ ] | `Manhattan()`; for probabilities equals `2 * total variation`. |
| Euclidean / L2 distance | [x] | [ ] | [x] | [ ] | [~] | [x] | [ ] | [x] | [ ] | [ ] | General distance frequently used on transformed data. |
| Chord distance | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | `Chord()` / `chord_distance`. |
| Hellinger distance | [x] | [ ] | [x] | [ ] | [~] | [~] | [x] | [ ] | [ ] | [ ] | `Hellinger()` / `hellinger_distance`; vegan and Microbiome.jl document related workflows. |
| Bhattacharyya coefficient / distance | [x] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [~] | [ ] | [ ] | `Bhattacharyya()`; probability-overlap family related to Hellinger. |
| Morisita / Morisita-Horn | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [~] | [~] | `MorisitaHorn()` / `morisita_horn_similarity`; vegan `morisita` and `horn`. |
| Canberra distance | [x] | [ ] | [x] | [ ] | [~] | [x] | [ ] | [x] | [ ] | [ ] | `Canberra()` / `canberra_distance`; denominator convention varies. |
| Clark / Cao / binomial / Gower / chi-square | [ ] | [ ] | [x] | [ ] | [~] | [~] | [ ] | [ ] | [ ] | [ ] | Ecological dissimilarities available elsewhere; not core here. |
| Kullback-Leibler divergence | [x] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [~] | [ ] | [x] | `KullbackLeibler()`; asymmetric with low-sample correction options. |
| Jeffreys divergence | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [~] | Symmetrized KL; not currently implemented. |
| Shannon entropy difference | [x] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [~] | `ShannonDifference()` compares entropy values, not shared taxa. |
| Jensen / Jensen-Shannon divergence / distance | [x] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [x] | [ ] | [~] | `JensenDifference()` and `JensenShannon()`; square root of JS divergence is a metric. |
| Aitchison / robust Aitchison | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [~] | Compositional-data distance after CLR-style transforms; requires zero handling. |
| Wasserstein / earth mover's | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [x] | [ ] | [ ] | Requires a ground distance among taxa/features. |
| UniFrac and phylogenetic pairwise distances | [ ] | [x] | [ ] | [ ] | [x] | [ ] | [~] | [ ] | [ ] | [~] | Out of scope here; Diversity.jl and scikit-bio are stronger phylogenetic options. |
| Multiple-site beta diversity | [ ] | [x] | [x] | [~] | [~] | [x] | [ ] | [ ] | [~] | [x] | Multi-community rather than pairwise-only comparison. |

## Reference

```@docs
Jaccard
SorensenDice
Overlap
BrayCurtis
Ruzicka
TotalVariation
Manhattan
Euclidean
Canberra
Hellinger
Chord
Bhattacharyya
KullbackLeibler
ShannonDifference
JensenDifference
JensenShannon
MorisitaHorn
similarity
dissimilarity
distance
labeled_similarity
labeled_dissimilarity
labeled_distance
jaccard_similarity
jaccard_index
jaccard_distance
sorensen_index
sorensen_dice_index
sorensen_distance
sorensen_dice_distance
sorensen_dice_dissimilarity
bray_curtis_distance
bray_curtis_dissimilarity
overlap_similarity
overlap_distance
ruzicka_similarity
quantitative_jaccard_similarity
ruzicka_distance
quantitative_jaccard_distance
total_variation_distance
manhattan_distance
euclidean_distance
canberra_distance
hellinger_distance
chord_distance
bhattacharyya_coefficient
bhattacharyya_distance
kullback_leibler_divergence
shannon_difference
jensen_difference
jensen_shannon_similarity
jensen_shannon_divergence
jensen_shannon_distance
morisita_horn_similarity
morisita_horn_distance
```
