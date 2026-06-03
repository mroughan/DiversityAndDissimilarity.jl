# Similarity And Dissimilarity Availability Checklist

This checklist compares pairwise similarity, dissimilarity, distance, and
overlap indices in `DiversityAndDissimilarity.jl` with a selection of commonly used
biodiversity, community-ecology, microbiome, and probability-comparison
packages. Alpha-diversity indices are covered separately in
[Diversity Index Availability Checklist](index-checklist.md).
For a broader domain catalog that includes indices not yet implemented here,
see [Similarity And Dissimilarity Catalog](similarity-index-catalog.md).

Last checked: 2026-05-15.

## Legend

- `[x]`: available directly or as a named/documented function or method.
- `[~]`: available indirectly, derivable from another output, or available with
  a different convention.
- `[ ]`: not documented as available in the checked source.

Package names and conventions are not perfectly standardized. Some packages
return a distance or dissimilarity where `DiversityAndDissimilarity.jl` exposes both
similarity and dissimilarity forms. This checklist intentionally focuses on
indices comparing sets/incidence data, abundance vectors, probability vectors,
or probability mass on a tree. It excludes edit distances and string metrics.

## Checklist

| Index or feature | DiversityAndDissimilarity.jl | Diversity.jl | vegan | iNEXT | scikit-bio | EcoPy | Microbiome.jl | SciPy | SpadeR | entropart / hill packages | Notes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Jaccard incidence similarity/distance | [x] | [ ] | [x] | [ ] | [~] | [x] | [x] | [x] | [x] | [ ] | `Jaccard()` / `jaccard_index`; vegan `vegdist(method="jaccard")`; EcoPy `distance(..., "jaccard")`; SciPy boolean `jaccard`; SpadeR includes incidence estimators. |
| Sorensen-Dice incidence similarity/distance | [x] | [ ] | [~] | [ ] | [~] | [x] | [ ] | [x] | [x] | [ ] | `SorensenDice()`; EcoPy `distance(..., "sorensen")`; SciPy boolean `dice`; vegan derives related binary forms through `betadiver` / transformations. |
| Overlap / Szymkiewicz-Simpson | [x] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [ ] | `Overlap()` / `overlap_similarity`; nestedness-sensitive incidence overlap. |
| Simple matching / Sokal-Michener | [ ] | [ ] | [ ] | [ ] | [~] | [x] | [ ] | [x] | [ ] | [ ] | Requires meaningful shared absences and a fixed species universe. |
| Russell-Rao | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [x] | [ ] | [ ] | Shared-absence-sensitive boolean coefficient. |
| Sokal-Sneath / Rogers-Tanimoto family | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [x] | [ ] | [ ] | Boolean coefficients available in SciPy; conventions vary. |
| Ochiai / binary cosine | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [~] | [ ] | [ ] | Cosine on binary incidence vectors; SciPy cosine can represent it after binary encoding. |
| Kulczynski incidence | [ ] | [ ] | [x] | [ ] | [~] | [ ] | [ ] | [ ] | [~] | [ ] | vegan includes Kulczynski among ecological dissimilarities; SpadeR includes related abundance/incidence estimators. |
| Mountford incidence | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | vegan `vegdist(method="mountford")`. |
| Raup-Crick incidence | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | Null-model incidence dissimilarity; depends on species pool/randomization assumptions. |
| Chao-Jaccard / Chao-Sorensen | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [x] | [~] | Adjusts overlap for unseen shared species; vegan has Chao dissimilarity; SpadeR focuses on estimated similarity/diversity. |
| Bray-Curtis dissimilarity | [x] | [ ] | [x] | [ ] | [~] | [x] | [x] | [x] | [~] | [ ] | `BrayCurtis()`; vegan `vegdist(method="bray")`; EcoPy `distance(..., "bray")`; Microbiome.jl `braycurtis`; SciPy `pdist(..., "braycurtis")`. |
| Quantitative Sorensen / Bray-Curtis similarity | [~] | [ ] | [~] | [ ] | [~] | [x] | [~] | [~] | [~] | [ ] | Derivable as `1 - bray_curtis_distance`; naming convention varies. |
| Ruzicka / quantitative Jaccard | [x] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [ ] | `Ruzicka()` / `ruzicka_similarity`; abundance version of Jaccard using `sum(min)/sum(max)`. |
| Percentage similarity / Renkonen | [ ] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | Probability-overlap index `sum(min(p,q))`; complement is total variation. |
| Total variation distance | [x] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [~] | [ ] | [ ] | `TotalVariation()` / `total_variation_distance`; equivalent to Bray-Curtis for normalized probability vectors. |
| Manhattan / L1 distance | [x] | [ ] | [x] | [ ] | [~] | [x] | [ ] | [x] | [ ] | [ ] | `Manhattan()` / `manhattan_distance`; for probabilities, equals `2 * total variation`. |
| Euclidean / L2 distance | [x] | [ ] | [x] | [ ] | [~] | [x] | [ ] | [x] | [ ] | [ ] | `Euclidean()` / `euclidean_distance`; general distance frequently used on transformed community/probability data. |
| Chord distance | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | `Chord()` / `chord_distance`; related to Euclidean distance after square-root normalization. |
| Hellinger distance | [x] | [ ] | [x] | [ ] | [~] | [~] | [x] | [ ] | [ ] | [ ] | `Hellinger()` / `hellinger_distance`; vegan `vegdist(method="hellinger")`; Microbiome.jl documents Hellinger beta diversity. |
| Bhattacharyya coefficient / distance | [x] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [~] | [ ] | [ ] | `Bhattacharyya()`; probability-overlap family related to Hellinger. |
| Morisita / Morisita-Horn | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [~] | [~] | `MorisitaHorn()` / `morisita_horn_similarity`; vegan `morisita` and `horn`; common abundance-overlap indices. |
| Horn information-theoretic overlap | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [~] | [x] | vegan `horn`; entropart/hill-family packages cover related entropy partitioning. |
| Kulczynski abundance | [ ] | [ ] | [x] | [ ] | [ ] | [x] | [ ] | [ ] | [~] | [ ] | vegan and EcoPy include Kulczynski-style ecological distances. |
| Canberra distance | [x] | [ ] | [x] | [ ] | [~] | [x] | [ ] | [x] | [ ] | [ ] | `Canberra()` / `canberra_distance`; relative abundance difference; denominator convention varies by implementation. |
| Clark distance | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | vegan `vegdist(method="clark")`. |
| Cao / CYd dissimilarity | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | vegan `vegdist(method="cao")`; designed for ecological count data. |
| Binomial deviance dissimilarity | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | vegan `vegdist(method="binomial")`. |
| Gower / altGower on abundance columns | [ ] | [ ] | [x] | [ ] | [~] | [x] | [ ] | [ ] | [ ] | [ ] | Include only when columns are species/probability variables, not arbitrary mixed metadata. |
| Chi-square distance | [ ] | [ ] | [x] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [ ] | vegan `vegdist(method="chisq")`; conventions vary. |
| Kullback-Leibler divergence | [x] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [~] | [ ] | [x] | `KullbackLeibler()` / `kullback_leibler_divergence`; asymmetric ``D_{KL}(left \Vert right)`` with Miller-Madow, pseudocount/shrinkage, and Good-Turing correction options. |
| Jeffreys divergence | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [~] | Symmetrized KL; not a metric. |
| Shannon entropy difference | [x] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [~] | `ShannonDifference()` / `shannon_difference`; compares entropy values rather than shared taxa. |
| Jensen / Jensen-Shannon divergence / distance | [x] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [x] | [ ] | [~] | `JensenDifference()` and `JensenShannon()`; includes low-sample correction options. SciPy includes Jensen-Shannon distance; square root of JS divergence is a metric. |
| Aitchison distance | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [~] | Compositional-data distance after CLR transform; requires zero handling. |
| Robust Aitchison distance | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [~] | Common in microbiome workflows for sparse compositions. |
| Earth mover's / Wasserstein on taxa probabilities | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [ ] | [x] | [ ] | [ ] | Requires a ground distance among taxa/features. |
| Unweighted UniFrac | [ ] | [x] | [ ] | [ ] | [x] | [ ] | [~] | [ ] | [ ] | [~] | scikit-bio directly documents unweighted UniFrac; Diversity.jl/entropart support phylogenetic diversity workflows. |
| Weighted UniFrac | [ ] | [x] | [ ] | [ ] | [x] | [ ] | [~] | [ ] | [ ] | [~] | Probability/abundance mass on phylogenetic branches. |
| Generalized UniFrac | [ ] | [ ] | [ ] | [ ] | [~] | [ ] | [~] | [ ] | [ ] | [ ] | Common microbiome distance family; implementation is usually in specialist packages. |
| Phylogenetic Jaccard / Sorensen | [ ] | [x] | [~] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [x] | Branch-set analogues of incidence overlap. |
| Jaccard/Sorensen turnover and nestedness partitioning | [ ] | [~] | [x] | [ ] | [ ] | [~] | [ ] | [ ] | [ ] | [~] | betapart is the specialist R package; vegan and some ecology packages expose related beta-diversity workflows. |
| Multiple-site beta diversity | [ ] | [x] | [x] | [~] | [~] | [x] | [ ] | [ ] | [~] | [x] | Multi-community rather than pairwise-only comparison. |

## Source Notes

- Diversity.jl: EcoJulia describes it as a package for measuring and
  partitioning diversity, including alpha, beta, and gamma diversity.
- vegan: `vegdist` documents Bray-Curtis, Jaccard, Hellinger, Euclidean,
  Manhattan, and many other ecological dissimilarities.
- iNEXT: focuses on interpolation/extrapolation and Hill-number diversity; it
  is not primarily a pairwise-distance package.
- scikit-bio: `skbio.diversity.beta_diversity` can use beta metrics, including
  SciPy-backed metrics, while `skbio.diversity.beta` directly documents UniFrac.
- EcoPy: `distance` documents Bray-Curtis, Jaccard, Sorensen, Euclidean,
  Manhattan, and other distances.
- Microbiome.jl: documents Bray-Curtis, Jaccard, and Hellinger beta diversity.
- SciPy: `scipy.spatial.distance.pdist` documents pairwise distances including
  Bray-Curtis, Jaccard, Dice, Euclidean, and cityblock/Manhattan.
- SpadeR: documents related similarity/dissimilarity estimation for abundance
  and incidence data.
- entropart: documents beta and gamma diversity, including phylogenetic and
  functional diversity, rather than a large catalogue of pairwise distances.

## References

- [EcoJulia Diversity.jl overview](https://libraries.io/julia/Diversity)
- [vegan vegdist](https://www.rdocumentation.org/packages/vegan/versions/1.8-6/topics/vegdist)
- [iNEXT package documentation](https://rdrr.io/cran/iNEXT/man/iNEXT-package.html)
- [scikit-bio community diversity](https://scikit.bio/docs/latest/diversity.html)
- [scikit-bio beta diversity](https://scikit.bio/docs/latest/generated/skbio.diversity.beta.html)
- [EcoPy matrix transformations and distance metrics](https://ecopy.readthedocs.io/en/latest/matrices.html)
- [Microbiome.jl diversity measures](https://docs.ecojulia.org/Microbiome.jl/v0.10/diversity/)
- [SciPy pairwise distances](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)
- [SpadeR package documentation](https://cran.r-universe.dev/SpadeR/doc/manual.html)
- [entropart package documentation](https://www.rdocumentation.org/packages/entropart/versions/1.6-16)
