# Diversity Index Availability Checklist

This checklist compares alpha-diversity indices and related single-assemblage
features in `DiversityIndices.jl` with a selection of commonly used
biodiversity and community ecology packages. Pairwise similarity and
dissimilarity indices are covered separately in
[Similarity And Dissimilarity Availability Checklist](similarity-checklist.md).

Last checked: 2026-05-13.

## Legend

- `[x]`: available directly or as a named/documented function or method.
- `[~]`: available indirectly, derivable from another output, or available with
  a different convention.
- `[ ]`: not documented as available in the checked source.

Package names and formulas are not perfectly standardized. In particular,
"Simpson" may mean ``\sum_i p_i^2``, ``1 - \sum_i p_i^2``, or
``1 / \sum_i p_i^2``, depending on the package. The notes call out these
convention differences where they matter.

## Checklist

| Index or feature | DiversityIndices.jl | Diversity.jl | vegan | iNEXT | scikit-bio | EcoPy | Microbiome.jl | SciPy | SpadeR | entropart | Notes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Observed richness / number of taxa | [x] | [x] | [x] | [x] | [x] | [x] | [ ] | [ ] | [x] | [x] | `Richness()` / `richness`; vegan `specnumber`; scikit-bio `sobs` / `observed_features`; EcoPy `spRich`. |
| Chao / ACE richness estimators | [x] | [ ] | [x] | [x] | [x] | [ ] | [ ] | [ ] | [x] | [~] | `Chao1()` / `chao1` and `ACE()` / `ace`; vegan `estimateR` / `specpool`; iNEXT `ChaoRichness`; scikit-bio `chao1` and `ace`; SpadeR is explicitly focused on richness prediction and diversity estimation. |
| Rarefied richness / rarefaction curves | [ ] | [ ] | [x] | [x] | [~] | [x] | [ ] | [ ] | [~] | [ ] | vegan `rarefy`; iNEXT interpolation/extrapolation curves; EcoPy `rarefy`. scikit-bio has related alpha tools but not the same ecology-facing rarefaction API in the checked page. |
| Shannon entropy, `H` | [x] | [~] | [x] | [~] | [x] | [x] | [x] | [ ] | [x] | [x] | `Shannon()` / `shannon_entropy`; defaults to `base=2` in this package. iNEXT focuses on Hill-number Shannon diversity; Diversity.jl emphasizes partitioned Hill/Jost measures. |
| Shannon effective diversity, ``\exp(H)`` or ``b^H`` | [x] | [x] | [~] | [x] | [x] | [x] | [ ] | [ ] | [x] | [x] | `effective_diversity(Shannon())`; vegan can derive it from `exp(diversity(..., "shannon"))`; scikit-bio `shannon(..., exp=true)`. |
| Shannon entropy estimation | [x] | [ ] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [x] | [x] | `Plugin()`, `MillerMadow()`, `HausserStrimmer()`, `Basharin()`, `AddGamma()`, and `ChaoShen()` are available here. iNEXT includes Chao-style Shannon estimation. entropart and SpadeR include estimation workflows. |
| Sample coverage | [x] | [ ] | [~] | [x] | [x] | [ ] | [ ] | [ ] | [x] | [x] | `SampleCoverage()` / `sample_coverage`; related coverage estimates are part of iNEXT, scikit-bio, SpadeR, and entropart workflows. |
| Simpson concentration, ``\sum_i p_i^2`` | [x] | [~] | [~] | [~] | [~] | [x] | [ ] | [ ] | [x] | [x] | This package's `Simpson()` returns concentration. vegan's `index="simpson"` returns ``1-D``; scikit-bio `simpson` also documents ``1-D``, while `dominance` / `simpson_d` cover dominance-style variants. |
| Gini-Simpson, ``1 - \sum_i p_i^2`` | [x] | [~] | [x] | [~] | [x] | [x] | [x] | [ ] | [x] | [~] | `GiniSimpson()`; vegan `diversity(..., "simpson")`; EcoPy `gini-simpson`; Microbiome.jl `ginisimpson`. |
| Greenberg / linguistic diversity index, ``1 - \sum_i p_i^2`` | [x] | [~] | [~] | [~] | [~] | [~] | [~] | [ ] | [~] | [~] | `GreenbergDiversityIndex()` / `LinguisticDiversityIndex()`; same formula as Gini-Simpson but interpreted as the probability that two people have different mother tongues. |
| Inverse Simpson, ``1 / \sum_i p_i^2`` | [x] | [x] | [x] | [x] | [x] | [x] | [ ] | [ ] | [x] | [x] | `InverseSimpson()`; vegan `index="invsimpson"`; scikit-bio `inv_simpson`; EcoPy number-equivalent output for Simpson. |
| Hill number, general order `q` | [x] | [x] | [~] | [x] | [x] | [~] | [ ] | [ ] | [x] | [x] | `Hill(q)`; Diversity.jl is centered on Hill/Jost diversity partitioning; vegan can use generalized profiles; iNEXT and scikit-bio expose Hill-number workflows. |
| Renyi entropy / diversity profile | [x] | [~] | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [~] | [x] | `Renyi(q; base=2)` / `renyi_entropy`; vegan `renyi`; scikit-bio `renyi`; entropart covers related generalized entropy/diversity profiles. |
| Tsallis entropy | [x] | [~] | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [x] | `Tsallis(q; base=2)` / `tsallis_entropy`; vegan `tsallis`; scikit-bio `tsallis`; entropart includes Tsallis-style entropy functions. |
| Fisher alpha | [x] | [ ] | [x] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [ ] | `FisherAlpha()` / `fisher_alpha`; vegan `fisher.alpha`; scikit-bio `fisher_alpha`. |
| Evenness / equitability | [x] | [ ] | [~] | [ ] | [x] | [x] | [ ] | [ ] | [ ] | [ ] | `PielouEvenness()` / `pielou_evenness`; vegan examples derive Pielou evenness from Shannon and richness; scikit-bio has Pielou, Simpson, Heip, and McIntosh evenness metrics; EcoPy `even`. |
| Dominance / Berger-Parker style metrics | [ ] | [ ] | [~] | [ ] | [x] | [x] | [ ] | [ ] | [ ] | [ ] | EcoPy `dominance`; scikit-bio lists dominance and Berger-Parker-related metrics. |
| UniFrac / Faith phylogenetic diversity | [ ] | [x] | [~] | [ ] | [x] | [ ] | [ ] | [ ] | [ ] | [x] | Diversity.jl supports phylogenetic diversity when used with Phylo.jl. scikit-bio has Faith PD and weighted/unweighted UniFrac. vegan has some tree-related tools but not the same UniFrac-centered API. |
| Rao quadratic entropy | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [ ] | [~] | entropart has phylogenetic entropy/diversity tooling; BAT and hilldiv are stronger specialist options for Rao-style measures. |
| Alpha/beta/gamma diversity partitioning | [ ] | [x] | [x] | [~] | [~] | [x] | [ ] | [ ] | [~] | [x] | Diversity.jl, vegan, EcoPy, and entropart provide broader partitioning workflows. SpadeR includes multi-community similarity/diversity estimation. |
| Functional diversity | [ ] | [~] | [~] | [ ] | [~] | [ ] | [ ] | [ ] | [ ] | [x] | entropart includes phylogenetic and functional diversity support; BAT, hillR, and hilldiv are also strong specialist options. |

## Source Notes

- Diversity.jl: EcoJulia describes it as a package for measuring and
  partitioning diversity, including alpha, beta, and gamma diversity; its
  phylogenetic page documents use with Phylo.jl for phylogenetic diversity.
- vegan: `diversity` documents Shannon, Simpson, inverse Simpson, Fisher alpha,
  and species richness; `specpool` / `estimateR` document extrapolated richness
  estimators including Chao and ACE.
- iNEXT: package documentation describes rarefaction/extrapolation for Hill
  numbers, especially species richness, Shannon diversity, and Simpson
  diversity, plus Chao-style estimators.
- scikit-bio: `skbio.diversity.alpha` lists richness, Shannon, Simpson,
  inverse Simpson, Hill, Renyi, Tsallis, evenness, dominance, and coverage
  metrics; `skbio.diversity.beta` directly documents UniFrac.
- EcoPy: `diversity` documents Shannon, Gini-Simpson, Simpson, dominance,
  species richness, and evenness; `rarefy` documents rarefied richness.
- Microbiome.jl: documents Shannon and Gini-Simpson alpha diversity plus
  Bray-Curtis, Jaccard, and Hellinger beta diversity.
- SciPy: `scipy.spatial.distance` is not an ecology-specific alpha-diversity
  package.
- SpadeR: documents species-richness prediction, diversity estimation, and
  related similarity/dissimilarity measures for abundance and incidence data.
- entropart: documents alpha, beta, and gamma diversity, including
  phylogenetic and functional diversity, with estimation-bias corrections.

## References

- [EcoJulia Diversity.jl overview](https://libraries.io/julia/Diversity)
- [Diversity.jl phylogenetic diversity](https://docs.ecojulia.org/Diversity.jl/dev/phylogenetics/)
- [vegan diversity](https://rdrr.io/cran/vegan/man/diversity.html)
- [vegan vegdist](https://www.rdocumentation.org/packages/vegan/versions/1.8-6/topics/vegdist)
- [vegan specpool / estimateR](https://rdrr.io/cran/vegan/man/specpool.html)
- [iNEXT package documentation](https://rdrr.io/cran/iNEXT/man/iNEXT-package.html)
- [scikit-bio community diversity](https://scikit.bio/docs/latest/diversity.html)
- [scikit-bio alpha diversity](https://scikit.bio/docs/latest/generated/skbio.diversity.alpha.html)
- [scikit-bio beta diversity](https://scikit.bio/docs/latest/generated/skbio.diversity.beta.html)
- [EcoPy species diversity](https://ecopy.readthedocs.io/en/latest/diversity.html)
- [EcoPy matrix transformations and distance metrics](https://ecopy.readthedocs.io/en/latest/matrices.html)
- [Microbiome.jl diversity measures](https://docs.ecojulia.org/Microbiome.jl/v0.10/diversity/)
- [SciPy pairwise distances](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html)
- [SpadeR package documentation](https://cran.r-universe.dev/SpadeR/doc/manual.html)
- [entropart package documentation](https://www.rdocumentation.org/packages/entropart/versions/1.6-16)
- [rdiversity package documentation](https://www.rdocumentation.org/packages/rdiversity/versions/1.2)
- [hillR documentation](https://phylotastic.r-universe.dev/hillR/doc/manual.html)
- [hilldiv package documentation](https://www.rdocumentation.org/packages/hilldiv/versions/1.5.1)
- [BAT package documentation](https://www.rdocumentation.org/packages/BAT/versions/2.11.1)
- [diverse package documentation](https://www.rdocumentation.org/packages/diverse/versions/0.1.1/topics/diversity)
