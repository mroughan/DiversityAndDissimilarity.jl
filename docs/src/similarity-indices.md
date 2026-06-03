# Similarity And Dissimilarity Indices

Similarity and dissimilarity indices compare two assemblages. In
`DiversityAndDissimilarity`, use [`similarity`](@ref), [`dissimilarity`](@ref), or the
alias [`distance`](@ref) with an index object:

```julia
similarity(Jaccard(), left, right)
dissimilarity(Jaccard(), left, right)
distance(Jaccard(), left, right)
similarity(SorensenDice(), left, right)
dissimilarity(BrayCurtis(), left, right)
distance(Hellinger(), left, right)
jensen_shannon_distance(left, right)
```

The two main styles are incidence comparisons and abundance comparisons.
See [Similarity And Dissimilarity Catalog](@ref) for a broader list of
set-based, abundance-based, and probability-based indices.

## Choosing A Measure

Use incidence indices when only presence/absence matters, or when abundance
counts are not comparable across samples. Use abundance indices when shared
species should count more if they are also similarly abundant. Use probability
distances when samples have different total effort and the comparison should be
about composition rather than absolute count.

Common defaults:

- Use [`Jaccard`](@ref) for a simple presence/absence turnover question.
- Use [`SorensenDice`](@ref) when shared presences should receive more weight
  than in Jaccard.
- Use [`Overlap`](@ref) when asking whether the smaller assemblage is nested in
  the larger one.
- Use [`BrayCurtis`](@ref) for ecological abundance data when abundance
  differences should matter and joint absences should not.
- Use [`Ruzicka`](@ref) when you want the abundance analogue of Jaccard.
- Use [`Hellinger`](@ref) or [`Chord`](@ref) before ordination-style workflows
  where Euclidean geometry on transformed relative abundances is desirable.
- Use [`TotalVariation`](@ref) for a direct probability-composition difference
  with a clear ``[0,1]`` range.
- Use [`JensenShannon`](@ref) for an information-theoretic, symmetric
  probability comparison; the square-root distance is metric.
- Use [`MorisitaHorn`](@ref) when dominant shared species should drive
  similarity and sample-size sensitivity should be reduced.

## Incidence Comparisons

Incidence comparisons use species presence or absence. [`Jaccard`](@ref),
[`SorensenDice`](@ref), and [`Overlap`](@ref) compare the species sets ``A`` and
``B``:

```math
J(A,B) = \\frac{|A \\cap B|}{|A \\cup B|}
```

```math
S(A,B) = \\frac{2|A \\cap B|}{|A| + |B|}
```

```math
O(A,B) = \\frac{|A \\cap B|}{\\min(|A|, |B|)}
```

| Index | Similarity range | Dissimilarity range | Most useful when |
|---|---:|---:|---|
| [`Jaccard`](@ref) | ``[0,1]`` | ``[0,1]`` | Measuring species turnover with no shared-absence term. |
| [`SorensenDice`](@ref) | ``[0,1]`` | ``[0,1]`` | Emphasizing shared species more strongly than Jaccard. |
| [`Overlap`](@ref) | ``[0,1]`` | ``[0,1]`` | Detecting nested assemblages, especially with unequal richness. |

For these indices, dissimilarity is defined as one minus similarity.

## Abundance Comparisons

Abundance comparisons use aligned abundance vectors ``x`` and ``y``. For
dictionaries, species are aligned by key. For numeric vectors with
`frequencies=true`, positions are treated as corresponding species and the
vectors must have the same length.

[`BrayCurtis`](@ref) and [`Ruzicka`](@ref) are both bounded ecological
dissimilarities that ignore joint absences:

```math
BC(x,y) = \\frac{\\sum_i |x_i - y_i|}{\\sum_i (x_i + y_i)}.
```

```math
R(x,y) = \\frac{\\sum_i \\min(x_i,y_i)}{\\sum_i \\max(x_i,y_i)}.
```

[`Canberra`](@ref) averages relative coordinate-wise differences:

```math
C(x,y) =
\\frac{1}{m}\\sum_{i:x_i+y_i>0}\\frac{|x_i-y_i|}{x_i+y_i}.
```

[`MorisitaHorn`](@ref) is an abundance-overlap similarity that is dominated by
common taxa:

```math
MH(x,y) =
\\frac{2\\sum_i x_i y_i}{(\\lambda_x + \\lambda_y) N_x N_y}.
```

| Index | Returned form | Range | Most useful when |
|---|---|---:|---|
| [`BrayCurtis`](@ref) | dissimilarity | ``[0,1]`` | Ecological abundance data with unequal sample totals; robust default for community matrices. |
| [`Ruzicka`](@ref) | similarity or ``1-R`` | ``[0,1]`` | You want a quantitative Jaccard interpretation using shared abundance over total abundance envelope. |
| [`Canberra`](@ref) | averaged distance | ``[0,1]`` | Rare taxa or low-abundance coordinates should have strong relative influence. |
| [`MorisitaHorn`](@ref) | similarity or ``1-MH`` | usually ``[0,1]`` | Dominant shared taxa should drive similarity; useful with count abundance data. |

## Probability Comparisons

Probability comparisons first normalize aligned abundance vectors to
probabilities ``p`` and ``q``. They are useful when sampling effort differs and
the scientific question is about composition.

```math
TV(p,q) = \\frac{1}{2}\\sum_i |p_i - q_i|
```

```math
H(p,q) =
\\frac{1}{\\sqrt{2}}\\sqrt{\\sum_i (\\sqrt{p_i} - \\sqrt{q_i})^2}
```

```math
JS(p,q) =
\\frac{1}{2}D_{KL}(p \\Vert m) + \\frac{1}{2}D_{KL}(q \\Vert m),
\\qquad m = \\frac{p+q}{2}.
```

```math
D_{KL}(p \\Vert q) = \\sum_i p_i \\log_b \\frac{p_i}{q_i}
```

```math
J_H(p,q) =
H_b\\left(\\frac{p+q}{2}\\right) - \\frac{H_b(p)+H_b(q)}{2}
```

| Index | Returned form | Range with default conventions | Most useful when |
|---|---|---:|---|
| [`TotalVariation`](@ref) | distance | ``[0,1]`` | You want the maximum compositional probability mass that differs between samples. |
| [`Manhattan`](@ref) | distance | ``[0,2]`` | You need an L1 distance; equals ``2TV`` for probability vectors. |
| [`Euclidean`](@ref) | distance | ``[0,\\sqrt{2}]`` | You need a familiar L2 geometry on relative abundances. |
| [`Hellinger`](@ref) | distance | ``[0,1]`` | You want a bounded metric that moderates dominant taxa through square-root transformation. |
| [`Chord`](@ref) | distance | ``[0,\\sqrt{2}]`` | You want Euclidean distance after square-root probability transformation. |
| [`Bhattacharyya`](@ref) | coefficient / distance | coefficient ``[0,1]``; distance ``[0,\\infty]`` | You want probability overlap, especially in classification or distributional overlap settings. |
| [`KullbackLeibler`](@ref) | asymmetric divergence | ``[0,\\infty]`` | You want ``D_{KL}(p \\Vert q)`` and have a directional reference/comparison interpretation. |
| [`ShannonDifference`](@ref) | dissimilarity | support-dependent | You want to compare entropy magnitudes; note this is **not** a distributional divergence — two assemblages with disjoint species but the same entropy score zero. |
| [`JensenDifference`](@ref) | divergence | ``[0,1]`` with `base=2` | You want the Jensen-Shannon divergence (raw, not square-rooted). Identical to `JensenShannon(distance=false)`. |
| [`JensenShannon`](@ref) | distance by default | ``[0,1]`` with `base=2` | You want a symmetric information-theoretic metric; the default square-root form is a proper distance. Use `distance=false` to get the raw divergence (same as `JensenDifference`). |

[`KullbackLeibler`](@ref) is directional: `dissimilarity(KullbackLeibler(), left, right)`
returns ``D_{KL}(left \\Vert right)``. It returns `Inf` if `right` has zero
probability where `left` has positive probability. Pairwise matrices for this
index are not symmetric: entry ``(i,j)`` is ``D_{KL}(i \\Vert j)`` and generally
differs from entry ``(j,i)``. You can check this programmatically with
`is_symmetric(KullbackLeibler())`, which returns `false`.

## Low-Sample Divergence Corrections

KL, Jensen difference, and Jensen-Shannon divergence/distance accept the same
estimator objects used for Shannon entropy:

```julia
kullback_leibler_divergence(left, right; estimator=MillerMadow())
kullback_leibler_divergence(left, right; estimator=AddGamma(1))    # Laplace
kullback_leibler_divergence(left, right; estimator=AddGamma(0.5))  # Jeffreys
kullback_leibler_divergence(left, right; estimator=HausserStrimmer())
kullback_leibler_divergence(left, right; estimator=ChaoShen())

jensen_difference(left, right; estimator=MillerMadow())
jensen_shannon_divergence(left, right; estimator=AddGamma(0.5), support=10)
jensen_shannon_distance(left, right; estimator=ChaoShen())
```

The correction modes have different interpretations:

- [`MillerMadow`](@ref) applies a first-order entropy bias correction: it
  subtracts ``(S-1)/(2n \log b)`` from the plugin divergence, where ``S`` and
  ``n`` are the observed support size and total count of the *left* distribution.
  This corrects for bias in ``H(p)`` but leaves the cross-entropy term
  uncorrected; it is a heuristic first-order adjustment.
- [`AddGamma`](@ref) applies Bayesian pseudocount smoothing. `AddGamma(1)` is
  Laplace smoothing and `AddGamma(0.5)` is Jeffreys smoothing.
- [`HausserStrimmer`](@ref) shrinks probabilities toward a uniform distribution
  over the aligned or supplied support.
- [`ChaoShen`](@ref) uses Good-Turing sample coverage to assign unseen mass to
  categories missing from one side and to a residual unseen category.

For `AddGamma` and `HausserStrimmer`, pass `support` when the finite category
universe is known and larger than the observed aligned support.

## Community Distance Matrices

Passing a community matrix as the only data argument computes all pairwise
comparisons across rows. Rows are sites/samples and columns are species/taxa.
The same method works for Tables.jl-compatible inputs, including DataFrames,
using `species` to choose the species columns.

```julia
distance(BrayCurtis(), community)
dissimilarity(Jaccard(), community)
similarity(SorensenDice(), community)
distance(Hellinger(), community)
jensen_shannon_distance(community)

bray_curtis_distance(table; species=[:oak, :ash, :elm])
```

Use labeled wrappers when sample/site identifiers should travel with the
matrix:

```julia
labeled_distance(BrayCurtis(), community; labels=["plot-a", "plot-b"])
labeled_distance(BrayCurtis(), table; label=:site, species=[:oak, :ash, :elm])
```

These return named tuples with `labels` and `matrix` fields.

Pairwise comparison of two vectors is linear in the aligned support size. A
full community distance matrix is more expensive: for ``M`` samples and ``P``
species columns, dense pairwise distances are ``O(M^2P)`` time and return an
``M \\times M`` result. This is often the limiting step for large datasets; see
[Scaling And Performance](@ref) for guidance on when to avoid building the full
matrix.

## Convenience Functions

Convenience functions call the same generic methods:

```julia
jaccard_index(left, right)
jaccard_similarity(left, right)
jaccard_distance(left, right)
sorensen_index(left, right)
sorensen_dice_index(left, right)
sorensen_distance(left, right)
sorensen_dice_distance(left, right)
sorensen_dice_dissimilarity(left, right)
bray_curtis_distance(left, right)
bray_curtis_dissimilarity(left, right)
overlap_similarity(left, right)
overlap_distance(left, right)
ruzicka_similarity(left, right)
quantitative_jaccard_similarity(left, right)
ruzicka_distance(left, right)
quantitative_jaccard_distance(left, right)
total_variation_distance(left, right)
manhattan_distance(left, right)
euclidean_distance(left, right)
canberra_distance(left, right)
hellinger_distance(left, right)
chord_distance(left, right)
bhattacharyya_coefficient(left, right)
bhattacharyya_distance(left, right)
kullback_leibler_divergence(left, right)
shannon_difference(left, right)
jensen_difference(left, right)
jensen_shannon_similarity(left, right)
jensen_shannon_divergence(left, right)
jensen_shannon_distance(left, right)
morisita_horn_similarity(left, right)
morisita_horn_distance(left, right)
```

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
