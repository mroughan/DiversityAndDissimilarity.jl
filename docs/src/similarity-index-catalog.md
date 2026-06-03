# Similarity And Dissimilarity Catalog

This catalog lists commonly used similarity, dissimilarity, distance, and
overlap indices whose inputs are either sets/incidence data or probability-like
abundance vectors. It intentionally excludes edit distances, string metrics, and
general-purpose distances whose primary interpretation is not ecological,
compositional, or probability based.

The package currently implements [`Jaccard`](@ref), [`SorensenDice`](@ref), and
[`BrayCurtis`](@ref). The remaining entries are useful candidates for future
implementation or for cross-package comparison.

## Notation

For incidence/set data, let ``A`` and ``B`` be the two observed species sets.
Equivalently, define:

- ``a = |A \cap B|``, shared presences.
- ``b = |A \setminus B|``, species only in ``A``.
- ``c = |B \setminus A|``, species only in ``B``.
- ``d``, shared absences, only meaningful when a fixed universe is defined.

For abundance or probability data, let ``x_i, y_i \ge 0`` be aligned abundance
vectors, and let ``p_i`` and ``q_i`` be the corresponding normalized
probabilities. Many quantitative indices are defined for abundances but become
probability distances after normalization.

## Incidence And Set Indices

| Index | Form | Similarity or dissimilarity | Notes |
|---|---|---:|---|
| Jaccard | ``a/(a+b+c)`` | similarity | Dissimilarity is ``(b+c)/(a+b+c)``. Also called Tanimoto for binary data in some fields. |
| Sorensen-Dice | ``2a/(2a+b+c)`` | similarity | Also called Dice, Sørensen, or Czekanowski for incidence data. |
| Overlap / Szymkiewicz-Simpson | ``a/\min(|A|, |B|)`` | similarity | Equals 1 when the smaller set is nested in the larger set. |
| Simpson turnover | ``\min(b,c)/(a+\min(b,c))`` | dissimilarity | Used in beta-diversity partitioning to isolate replacement/turnover. |
| Nestedness-resultant Jaccard | ``\beta_{jac} - \beta_{jtu}`` | dissimilarity component | Jaccard-family nestedness component used by betapart-style decompositions. |
| Nestedness-resultant Sorensen | ``\beta_{sor} - \beta_{sim}`` | dissimilarity component | Sorensen-family nestedness component used by betapart-style decompositions. |
| Simple matching | ``(a+d)/(a+b+c+d)`` | similarity | Uses shared absences, so it depends on a meaningful universe. |
| Russell-Rao | ``a/(a+b+c+d)`` | similarity | Also depends on a fixed universe. |
| Sokal-Sneath | ``a/(a+2(b+c))`` | similarity | Several Sokal-Sneath conventions exist; document the convention used. |
| Rogers-Tanimoto | ``(a+d)/(a+d+2(b+c))`` | similarity | Shared absences matter. |
| Ochiai / binary cosine | ``a/\sqrt{(a+b)(a+c)}`` | similarity | Cosine similarity for binary incidence vectors. |
| Braun-Blanquet | ``a/\max(|A|, |B|)`` | similarity | Penalizes unmatched richness by the larger set size. |
| Kulczynski incidence | ``\frac{1}{2}\left(a/(a+b) + a/(a+c)\right)`` | similarity | Average of conditional overlap in each direction. |
| Mountford | implicit function of ``a,b,c`` | dissimilarity or similarity | Incidence index available in vegan; conventions are less transparent than Jaccard/Sorensen. |
| Raup-Crick | probability from null co-occurrence model | dissimilarity | Incidence-based probabilistic index; depends on the null model and species pool. |
| Chao-Jaccard | abundance-adjusted incidence overlap | similarity | Adjusts Jaccard-style overlap for unseen shared species. |
| Chao-Sorensen | abundance-adjusted incidence overlap | similarity | Adjusts Sorensen-style overlap for unseen shared species. |

## Abundance And Probability Overlap

| Index | Form | Similarity or dissimilarity | Notes |
|---|---|---:|---|
| Bray-Curtis | ``\sum_i |x_i-y_i| / \sum_i(x_i+y_i)`` | dissimilarity | For normalized probabilities, equals total variation distance. |
| Quantitative Sorensen | ``2\sum_i \min(x_i,y_i)/(\sum_i x_i+\sum_i y_i)`` | similarity | Complement of Bray-Curtis under the usual abundance convention. |
| Ruzicka / quantitative Jaccard | ``\sum_i \min(x_i,y_i)/\sum_i \max(x_i,y_i)`` | similarity | Also called abundance Jaccard. |
| Percentage similarity / Renkonen | ``\sum_i \min(p_i,q_i)`` | similarity | Complement is total variation distance. |
| Total variation | ``\frac{1}{2}\sum_i |p_i-q_i|`` | distance | Probability distance; equivalent to Bray-Curtis for normalized vectors. |
| Manhattan / L1 | ``\sum_i |p_i-q_i|`` | distance | Twice total variation for probabilities. |
| Euclidean / L2 | ``\sqrt{\sum_i(p_i-q_i)^2}`` | distance | Useful for probabilities but not compositionally invariant. |
| Chord | ``\sqrt{\sum_i(\sqrt{p_i}-\sqrt{q_i})^2}`` | distance | Closely related to Hellinger distance; common after square-root transformation. |
| Hellinger | ``\frac{1}{\sqrt{2}}\sqrt{\sum_i(\sqrt{p_i}-\sqrt{q_i})^2}`` | distance | Metric on probability distributions. |
| Bhattacharyya coefficient | ``\sum_i\sqrt{p_iq_i}`` | similarity | Hellinger and Bhattacharyya distances are transformations of this coefficient. |
| Bhattacharyya distance | ``-\log\sum_i\sqrt{p_iq_i}`` | distance-like divergence | Not always finite when supports are disjoint. |
| Morisita-Horn | overlap based on ``\sum_i x_i y_i`` and Simpson concentration | similarity | Abundance overlap that downweights richness differences relative to dominant taxa. |
| Horn | information-theoretic overlap | similarity/dissimilarity | Used in vegan and ecological overlap literature. |
| Kulczynski abundance | average of shared abundance fractions | similarity/dissimilarity | Vegan exposes a quantitative Kulczynski dissimilarity. |
| Canberra | average of ``|x_i-y_i|/(x_i+y_i)`` terms | distance | Strongly weights rare species; denominator convention varies by package. |
| Clark | Euclidean form of Canberra-style relative differences | distance | Available in vegan. |
| Cao / CYd | log-transformed abundance dissimilarity | dissimilarity | Designed for ecological count data with many zeros. |
| Binomial deviance | binomial likelihood-style abundance dissimilarity | dissimilarity | Available in vegan. |
| Gower abundance | mean absolute difference after scaling | dissimilarity | Include only when variables are species abundance/probability columns. |

## Probability Divergences

These compare normalized probability vectors. Some are true metrics, while
others are asymmetric divergences.

| Index | Form | Type | Notes |
|---|---|---:|---|
| Kullback-Leibler | ``\sum_i p_i\log(p_i/q_i)`` | divergence | Asymmetric; returns infinity when ``p_i > 0`` and ``q_i = 0``. |
| Jeffreys | ``D_{KL}(p\|q)+D_{KL}(q\|p)`` | divergence | Symmetric but not a metric. |
| Shannon difference | ``|H(p)-H(q)|`` | dissimilarity | Compares entropy values rather than species overlap. |
| Jensen difference | ``H((p+q)/2) - (H(p)+H(q))/2`` | divergence | For Shannon entropy this is Jensen-Shannon divergence. |
| Jensen-Shannon | ``\frac{1}{2}D_{KL}(p\|m)+\frac{1}{2}D_{KL}(q\|m)``, ``m=(p+q)/2`` | divergence | Square root is a metric. |
| Jensen-Shannon distance | ``\sqrt{JS(p,q)}`` | distance | Bounded and symmetric. |
| Chi-square | ``\sum_i (p_i-q_i)^2/(p_i+q_i)`` or related conventions | distance/divergence | Multiple ecological and statistical conventions exist. |
| Hellinger | see above | distance | Often preferred over raw Euclidean distance for community composition. |
| Aitchison | Euclidean distance after centered log-ratio transform | distance | For compositional data; needs zero handling. |
| Robust Aitchison | robust CLR-based distance | distance | Common in microbiome workflows with sparse compositions. |
| Earth mover's / Wasserstein | minimum transport cost between probability masses | distance | Requires a ground distance among taxa/features. |

## Phylogenetic Set And Probability Distances

These compare sets or probabilities distributed on a tree rather than only on
taxon labels.

| Index | Input | Type | Notes |
|---|---|---:|---|
| Unweighted UniFrac | incidence on phylogenetic branches | distance | Branch-length-weighted set dissimilarity. |
| Weighted UniFrac | abundance/probability mass on branches | distance | Incorporates relative abundance along the tree. |
| Generalized UniFrac | abundance/probability mass on branches | distance family | Tuning parameter controls the influence of abundant versus rare lineages. |
| Phylogenetic Jaccard | branch sets | dissimilarity | Jaccard-style comparison over phylogenetic branch sets. |
| Phylogenetic Sorensen | branch sets | dissimilarity | Sorensen-style comparison over phylogenetic branch sets. |

## Practical Selection Notes

- Use incidence/set indices when only presence/absence is meaningful, or when
  sampling effort makes abundance comparisons unsafe.
- Use Jaccard or Sorensen when shared absences are not meaningful. Use simple
  matching or Rogers-Tanimoto only with a well-defined species universe.
- Use Bray-Curtis, quantitative Sorensen, Ruzicka, or Renkonen for ecological
  abundance vectors where total abundance and shared composition matter.
- Use Hellinger/chord-style distances when Euclidean methods are desired but raw
  abundance differences overemphasize dominant species.
- Use Jensen-Shannon, Hellinger, total variation, or Aitchison-style distances
  when the object of comparison is explicitly a probability composition.
- Use UniFrac-style distances only when a phylogeny or feature tree is part of
  the model.
- Always document whether a function returns similarity, dissimilarity, or
  distance, because many names are used for both a quantity and its complement.

## Sources

- [vegan `vegdist`](https://vegandevs.github.io/vegan/reference/vegdist.html)
- [vegan `vegdist` source methods](https://rdrr.io/cran/vegan/src/R/vegdist.R)
- [scikit-bio community diversity](https://scikit.bio/docs/latest/diversity.html)
- [scikit-bio beta diversity measures](https://scikit.bio/docs/dev/generated/skbio.diversity.beta.html)
- [betapart incidence beta diversity](https://search.r-project.org/CRAN/refmans/betapart/html/beta.pair.html)
- [betapart package overview](https://andres-baselga.r-universe.dev/betapart)
- [SpadeR package manual](https://cran.r-universe.dev/SpadeR/doc/manual.html)
